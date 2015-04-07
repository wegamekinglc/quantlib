/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Klaus Spanderen
 Copyright (C) 2015 Peter Caspers

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file levenbergmarquardt.hpp
    \brief Levenberg-Marquardt optimization method
*/

#ifndef quantlib_optimization_levenberg_marquardt_hpp
#define quantlib_optimization_levenberg_marquardt_hpp

#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/lmdif.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif
#include <boost/bind.hpp>
#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic pop
#endif

namespace QuantLib {

//! Levenberg-Marquardt optimization method
/*! This implementation is based on MINPACK
    (<http://www.netlib.org/minpack>,
    <http://www.netlib.org/cephes/linalg.tgz>)
*/
template <class T> class LevenbergMarquardt_t : public OptimizationMethod_t<T> {
  public:
    LevenbergMarquardt_t(T epsfcn = 1.0e-8, T xtol = 1.0e-8, T gtol = 1.0e-8);
    virtual typename EndCriteria_t<T>::Type
    minimize(Problem_t<T> &P,
             const EndCriteria_t<T> &endCriteria //= EndCriteria_t<T>()
             );
    //      = EndCriteria_t<T>(400, 1.0e-8, 1.0e-8)
    virtual Integer getInfo() const;
    void fcn(int m, int n, T *x, T *fvec, int *iflag);

  private:
    Problem_t<T> *currentProblem_;
    Array_t<T> initCostValues_;
    mutable Integer info_;
    const T epsfcn_, xtol_, gtol_;
};

typedef LevenbergMarquardt_t<Real> LevenbergMarquardt;

// implementation

template <class T>
LevenbergMarquardt_t<T>::LevenbergMarquardt_t(T epsfcn, T xtol, T gtol)
    : info_(0), epsfcn_(epsfcn), xtol_(xtol), gtol_(gtol) {}

template <class T> Integer LevenbergMarquardt_t<T>::getInfo() const {
    return info_;
}

template <class T>
typename EndCriteria_t<T>::Type
LevenbergMarquardt_t<T>::minimize(Problem_t<T> &P,
                                  const EndCriteria_t<T> &endCriteria) {
    EndCriteria::Type ecType = EndCriteria::None;
    P.reset();
    Array_t<T> x_ = P.currentValue();
    currentProblem_ = &P;
    initCostValues_ = P.costFunction().values(x_);
    int m = initCostValues_.size();
    int n = x_.size();
    boost::scoped_array<T> xx(new double[n]);
    std::copy(x_.begin(), x_.end(), xx.get());
    boost::scoped_array<T> fvec(new double[m]);
    boost::scoped_array<T> diag(new double[n]);
    int mode = 1;
    T factor = 1;
    int nprint = 0;
    int info = 0;
    int nfev = 0;
    boost::scoped_array<T> fjac(new double[m * n]);
    int ldfjac = m;
    boost::scoped_array<int> ipvt(new int[n]);
    boost::scoped_array<T> qtf(new double[n]);
    boost::scoped_array<T> wa1(new double[n]);
    boost::scoped_array<T> wa2(new double[n]);
    boost::scoped_array<T> wa3(new double[n]);
    boost::scoped_array<T> wa4(new double[m]);
    // requirements; check here to get more detailed error messages.
    QL_REQUIRE(n > 0, "no variables given");
    QL_REQUIRE(m >= n, "less functions (" << m << ") than available variables ("
                                          << n << ")");
    QL_REQUIRE(endCriteria.functionEpsilon() >= 0.0, "negative f tolerance");
    QL_REQUIRE(xtol_ >= 0.0, "negative x tolerance");
    QL_REQUIRE(gtol_ >= 0.0, "negative g tolerance");
    QL_REQUIRE(endCriteria.maxIterations() > 0, "null number of evaluations");

    // call lmdif to minimize the sum of the squares of m functions
    // in n variables by the Levenberg-Marquardt algorithm.
    typename MINPACK::LmdifCostFunction_t<T>::Type lmdifCostFunction =
        boost::bind(&LevenbergMarquardt_t<T>::fcn, this, _1, _2, _3, _4, _5);
    MINPACK::lmdif(
        m, n, xx.get(), fvec.get(),
        static_cast<T>(endCriteria.functionEpsilon()), static_cast<T>(xtol_),
        static_cast<T>(gtol_), static_cast<int>(endCriteria.maxIterations()),
        static_cast<T>(epsfcn_), diag.get(), mode, factor, nprint, &info, &nfev,
        fjac.get(), ldfjac, ipvt.get(), qtf.get(), wa1.get(), wa2.get(),
        wa3.get(), wa4.get(), lmdifCostFunction);
    info_ = info;
    // check requirements & endCriteria evaluation
    QL_REQUIRE(info != 0, "MINPACK: improper input parameters");
    // QL_REQUIRE(info != 6, "MINPACK: ftol is too small. no further "
    //                               "reduction in the sum of squares "
    //                               "is possible.");
    if (info != 6)
        ecType = QuantLib::EndCriteria::StationaryFunctionValue;
    // QL_REQUIRE(info != 5, "MINPACK: number of calls to fcn has "
    //                               "reached or exceeded maxfev.");
    endCriteria.checkMaxIterations(nfev, ecType);
    QL_REQUIRE(info != 7, "MINPACK: xtol is too small. no further "
                          "improvement in the approximate "
                          "solution x is possible.");
    QL_REQUIRE(info != 8, "MINPACK: gtol is too small. fvec is "
                          "orthogonal to the columns of the "
                          "jacobian to machine precision.");
    // set problem
    std::copy(xx.get(), xx.get() + n, x_.begin());
    P.setCurrentValue(x_);
    P.setFunctionValue(P.costFunction().value(x_));

    return ecType;
}

template <class T>
void LevenbergMarquardt_t<T>::fcn(int, int n, T *x, T *fvec, int *) {
    Array_t<T> xt(n);
    std::copy(x, x + n, xt.begin());
    // constraint handling needs some improvement in the future:
    // starting point should not be close to a constraint violation
    if (currentProblem_->constraint().test(xt)) {
        const Array_t<T> &tmp = currentProblem_->values(xt);
        std::copy(tmp.begin(), tmp.end(), fvec);
    } else {
        std::copy(initCostValues_.begin(), initCostValues_.end(), fvec);
    }
}

} // namespace QuantLib

#endif
