/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Nicolas Di Césaré
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

/*! \file costfunction.hpp
    \brief Optimization cost function class
*/

#ifndef quantlib_optimization_costfunction_h
#define quantlib_optimization_costfunction_h

#include <ql/math/array.hpp>

namespace QuantLib {

//!  Cost function abstract class for optimization problem
template <class T> class CostFunction_t {
  public:
    virtual ~CostFunction_t() {}
    //! method to overload to compute the cost function value in x
    virtual T value(const Array_t<T> &x) const = 0;
    //! method to overload to compute the cost function values in x
    virtual Disposable<Array_t<T> > values(const Array_t<T> &x) const = 0;

    //! method to overload to compute grad_f, the first derivative of
    //  the cost function with respect to x
    virtual void gradient(Array_t<T> &grad, const Array_t<T> &x) const {
        T eps = finiteDifferenceEpsilon(), fp, fm;
        Array_t<T> xx(x);
        for (Size i = 0; i < x.size(); i++) {
            xx[i] += eps;
            fp = value(xx);
            xx[i] -= 2.0 * eps;
            fm = value(xx);
            grad[i] = 0.5 * (fp - fm) / eps;
            xx[i] = x[i];
        }
    }

    //! method to overload to compute grad_f, the first derivative of
    //  the cost function with respect to x and also the cost function
    virtual T valueAndGradient(Array_t<T> &grad, const Array_t<T> &x) const {
        gradient(grad, x);
        return value(x);
    }

    //! Default epsilon for finite difference method :
    virtual T finiteDifferenceEpsilon() const { return 1e-8; }
};

typedef CostFunction_t<Real> CostFunction;

template <class T> class ParametersTransformation_t {
  public:
    virtual ~ParametersTransformation_t() {}
    virtual Array_t<T> direct(const Array_t<T> &x) const = 0;
    virtual Array_t<T> inverse(const Array_t<T> &x) const = 0;
};

typedef ParametersTransformation_t<Real> ParametersTransformation;
}

#endif
