/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl
 Copyright (C) 2011 Ferdinando Ametrano
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

/*! \file tridiagonaloperator.hpp
    \brief tridiagonal operator
*/

#ifndef quantlib_tridiagonal_operator_hpp
#define quantlib_tridiagonal_operator_hpp

#include <ql/math/array.hpp>
#include <ql/math/comparison.hpp>
#include <boost/shared_ptr.hpp>

namespace QuantLib {

//! Base implementation for tridiagonal operator
/*! \warning to use real time-dependant algebra, you must overload
             the corresponding operators in the inheriting
             time-dependent class.

    \ingroup findiff
*/
template <class T> class TridiagonalOperator_t {
    // unary operators
    template<class S> friend Disposable<TridiagonalOperator_t<S> > operator+
        (const TridiagonalOperator_t<S> &);
    template<class S> friend Disposable<TridiagonalOperator_t<S> > operator-
        (const TridiagonalOperator_t<S> &);
    // binary operators
    template<class S> friend Disposable<TridiagonalOperator_t<S> > operator+
        (const TridiagonalOperator_t<S> &, const TridiagonalOperator_t<S> &);
    template<class S> friend Disposable<TridiagonalOperator_t<S> > operator-
        (const TridiagonalOperator_t<S> &, const TridiagonalOperator_t<S> &);
    template<class S> friend Disposable<TridiagonalOperator_t<S> >
    operator*(T, const TridiagonalOperator_t<S> &);
    template<class S> friend Disposable<TridiagonalOperator_t<S> >
    operator*(const TridiagonalOperator_t<S> &, S);
    template<class S> friend Disposable<TridiagonalOperator_t<S> > operator/
        (const TridiagonalOperator_t<S> &, T);

  public:
    typedef Array_t<T> array_type;
    // constructors
    explicit TridiagonalOperator_t(Size size = 0);
    TridiagonalOperator_t(const Array_t<T> &low, const Array_t<T> &mid,
                          const Array_t<T> &high);
    TridiagonalOperator_t(const Disposable<TridiagonalOperator_t<T> > &);
    TridiagonalOperator_t &
    operator=(const Disposable<TridiagonalOperator_t<T> > &);
    //! \name Operator interface
    //@{
    //! apply operator to a given array
    Disposable<Array_t<T> > applyTo(const Array_t<T> &v) const;
    //! solve linear system for a given right-hand side
    Disposable<Array_t<T> > solveFor(const Array_t<T> &rhs) const;
    /*! solve linear system for a given right-hand side
        without result Array_t<T> allocation. The rhs and result parameters
        can be the same Array, in which case rhs will be changed
    */
    void solveFor(const Array_t<T> &rhs, Array_t<T> &result) const;
    //! solve linear system with SOR approach
    Disposable<Array_t<T> > SOR(const Array_t<T> &rhs, T tol) const;
    //! identity instance
    static Disposable<TridiagonalOperator_t<T> > identity(Size size);
    //@}
    //! \name Inspectors
    //@{
    Size size() const { return n_; }
    bool isTimeDependent() const { return !!timeSetter_; }
    const Array_t<T> &lowerDiagonal() const { return lowerDiagonal_; }
    const Array_t<T> &diagonal() const { return diagonal_; }
    const Array_t<T> &upperDiagonal() const { return upperDiagonal_; }
    //@}
    //! \name Modifiers
    //@{
    void setFirstRow(T, T);
    void setMidRow(Size, T, T, T);
    void setMidRows(T, T, T);
    void setLastRow(T, T);
    void setTime(Time t);
    //@}
    //! \name Utilities
    //@{
    void swap(TridiagonalOperator_t<T> &);
    //@}
    //! encapsulation of time-setting logic
    class TimeSetter {
      public:
        virtual ~TimeSetter() {}
        virtual void setTime(Time t, TridiagonalOperator_t<T> &L) const = 0;
    };

  protected:
    Size n_;
    Array_t<T> diagonal_, lowerDiagonal_, upperDiagonal_;
    mutable Array_t<T> temp_;
    boost::shared_ptr<TimeSetter> timeSetter_;
};

typedef TridiagonalOperator_t<Real> TridiagonalOperator;

/* \relates TridiagonalOperator_t */
template <class T>
void swap(TridiagonalOperator_t<T> &, TridiagonalOperator_t<T> &);

// inline definitions

template <class T>
inline TridiagonalOperator_t<T> &TridiagonalOperator_t<T>::
operator=(const Disposable<TridiagonalOperator_t<T> > &from) {
    swap(const_cast<Disposable<TridiagonalOperator_t<T> > &>(from));
    return *this;
}

template <class T>
inline void TridiagonalOperator_t<T>::setFirstRow(T valB, T valC) {
    diagonal_[0] = valB;
    upperDiagonal_[0] = valC;
}

template <class T>
inline void TridiagonalOperator_t<T>::setMidRow(Size i, T valA, T valB,
                                                T valC) {
    QL_REQUIRE(i >= 1 && i <= n_ - 2,
               "out of range in TridiagonalSystem::setMidRow");
    lowerDiagonal_[i - 1] = valA;
    diagonal_[i] = valB;
    upperDiagonal_[i] = valC;
}

template <class T>
inline void TridiagonalOperator_t<T>::setMidRows(T valA, T valB, T valC) {
    for (Size i = 1; i <= n_ - 2; i++) {
        lowerDiagonal_[i - 1] = valA;
        diagonal_[i] = valB;
        upperDiagonal_[i] = valC;
    }
}

template <class T>
inline void TridiagonalOperator_t<T>::setLastRow(T valA, T valB) {
    lowerDiagonal_[n_ - 2] = valA;
    diagonal_[n_ - 1] = valB;
}

template <class T> inline void TridiagonalOperator_t<T>::setTime(Time t) {
    if (timeSetter_)
        timeSetter_->setTime(t, *this);
}

template <class T>
inline void TridiagonalOperator_t<T>::swap(TridiagonalOperator_t<T> &from) {
    using std::swap;
    swap(n_, from.n_);
    diagonal_.swap(from.diagonal_);
    lowerDiagonal_.swap(from.lowerDiagonal_);
    upperDiagonal_.swap(from.upperDiagonal_);
    temp_.swap(from.temp_);
    swap(timeSetter_, from.timeSetter_);
}

// Time constant algebra

template <class T>
inline Disposable<TridiagonalOperator_t<T> >
operator+(const TridiagonalOperator_t<T> &D) {
    TridiagonalOperator_t<T> D1 = D;
    return D1;
}

template <class T>
inline Disposable<TridiagonalOperator_t<T> >
operator-(const TridiagonalOperator_t<T> &D) {
    Array_t<T> low = -D.lowerDiagonal_, mid = -D.diagonal_,
               high = -D.upperDiagonal_;
    TridiagonalOperator_t<T> result(low, mid, high);
    return result;
}

template <class T>
inline Disposable<TridiagonalOperator_t<T> >
operator+(const TridiagonalOperator_t<T> &D1,
          const TridiagonalOperator_t<T> &D2) {
    Array_t<T> low = D1.lowerDiagonal_ + D2.lowerDiagonal_,
               mid = D1.diagonal_ + D2.diagonal_,
               high = D1.upperDiagonal_ + D2.upperDiagonal_;
    TridiagonalOperator_t<T> result(low, mid, high);
    return result;
}

template <class T>
inline Disposable<TridiagonalOperator_t<T> >
operator-(const TridiagonalOperator_t<T> &D1,
          const TridiagonalOperator_t<T> &D2) {
    Array_t<T> low = D1.lowerDiagonal_ - D2.lowerDiagonal_,
               mid = D1.diagonal_ - D2.diagonal_,
               high = D1.upperDiagonal_ - D2.upperDiagonal_;
    TridiagonalOperator_t<T> result(low, mid, high);
    return result;
}

template <class T>
inline Disposable<TridiagonalOperator_t<T> >
operator*(T a, const TridiagonalOperator_t<T> &D) {
    Array_t<T> low = D.lowerDiagonal_ * a, mid = D.diagonal_ * a,
               high = D.upperDiagonal_ * a;
    TridiagonalOperator_t<T> result(low, mid, high);
    return result;
}

template <class T>
inline Disposable<TridiagonalOperator_t<T> >
operator*(const TridiagonalOperator_t<T> &D, T a) {
    Array_t<T> low = D.lowerDiagonal_ * a, mid = D.diagonal_ * a,
               high = D.upperDiagonal_ * a;
    TridiagonalOperator_t<T> result(low, mid, high);
    return result;
}

template <class T>
inline Disposable<TridiagonalOperator_t<T> >
operator/(const TridiagonalOperator_t<T> &D, T a) {
    Array_t<T> low = D.lowerDiagonal_ / a, mid = D.diagonal_ / a,
               high = D.upperDiagonal_ / a;
    TridiagonalOperator_t<T> result(low, mid, high);
    return result;
}

template <class T>
inline void swap(TridiagonalOperator_t<T> &L1, TridiagonalOperator_t<T> &L2) {
    L1.swap(L2);
}

// implementation

template <class T> TridiagonalOperator_t<T>::TridiagonalOperator_t(Size size) {
    if (size >= 2) {
        n_ = size;
        diagonal_ = Array_t<T>(size);
        lowerDiagonal_ = Array_t<T>(size - 1);
        upperDiagonal_ = Array_t<T>(size - 1);
        temp_ = Array_t<T>(size);
    } else if (size == 0) {
        n_ = 0;
        diagonal_ = Array_t<T>(0);
        lowerDiagonal_ = Array_t<T>(0);
        upperDiagonal_ = Array_t<T>(0);
        temp_ = Array_t<T>(0);
    } else {
        QL_FAIL("invalid size (" << size << ") for tridiagonal operator "
                                            "(must be null or >= 2)");
    }
}

template <class T>
TridiagonalOperator_t<T>::TridiagonalOperator_t(const Array_t<T> &low,
                                                const Array_t<T> &mid,
                                                const Array_t<T> &high)
    : n_(mid.size()), diagonal_(mid), lowerDiagonal_(low), upperDiagonal_(high),
      temp_(n_) {
    QL_REQUIRE(low.size() == n_ - 1, "low diagonal vector of size "
                                         << low.size() << " instead of "
                                         << n_ - 1);
    QL_REQUIRE(high.size() == n_ - 1, "high diagonal vector of size "
                                          << high.size() << " instead of "
                                          << n_ - 1);
}

template <class T>
TridiagonalOperator_t<T>::TridiagonalOperator_t(
    const Disposable<TridiagonalOperator_t<T> > &from) {
    swap(const_cast<Disposable<TridiagonalOperator_t<T> > &>(from));
}

template <class T>
Disposable<Array_t<T> >
TridiagonalOperator_t<T>::applyTo(const Array_t<T> &v) const {
    QL_REQUIRE(n_ != 0, "uninitialized TridiagonalOperator");
    QL_REQUIRE(v.size() == n_, "vector of the wrong size "
                                   << v.size() << " instead of " << n_);
    Array_t<T> result(n_);
    std::transform(diagonal_.begin(), diagonal_.end(), v.begin(),
                   result.begin(), std::multiplies<T>());

    // matricial product
    result[0] += upperDiagonal_[0] * v[1];
    for (Size j = 1; j <= n_ - 2; j++)
        result[j] +=
            lowerDiagonal_[j - 1] * v[j - 1] + upperDiagonal_[j] * v[j + 1];
    result[n_ - 1] += lowerDiagonal_[n_ - 2] * v[n_ - 2];

    return result;
}

template <class T>
Disposable<Array_t<T> >
TridiagonalOperator_t<T>::solveFor(const Array_t<T> &rhs) const {

    Array_t<T> result(rhs.size());
    solveFor(rhs, result);
    return result;
}

template <class T>
void TridiagonalOperator_t<T>::solveFor(const Array_t<T> &rhs,
                                        Array_t<T> &result) const {

    QL_REQUIRE(n_ != 0, "uninitialized TridiagonalOperator");
    QL_REQUIRE(rhs.size() == n_, "rhs vector of size " << rhs.size()
                                                       << " instead of " << n_);

    T bet = diagonal_[0];
    QL_REQUIRE(!close(bet, T(0.0)), "diagonal's first element ("
                                     << bet << ") cannot be close to zero");
    result[0] = rhs[0] / bet;
    for (Size j = 1; j <= n_ - 1; ++j) {
        temp_[j] = upperDiagonal_[j - 1] / bet;
        bet = diagonal_[j] - lowerDiagonal_[j - 1] * temp_[j];
        QL_ENSURE(!close(bet, T(0.0)), "division by zero");
        result[j] = (rhs[j] - lowerDiagonal_[j - 1] * result[j - 1]) / bet;
    }
    // cannot be j>=0 with Size j
    for (Size j = n_ - 2; j > 0; --j)
        result[j] -= temp_[j + 1] * result[j + 1];
    result[0] -= temp_[1] * result[1];
}

template <class T>
Disposable<Array_t<T> > TridiagonalOperator_t<T>::SOR(const Array_t<T> &rhs,
                                                      T tol) const {
    QL_REQUIRE(n_ != 0, "uninitialized TridiagonalOperator");
    QL_REQUIRE(rhs.size() == n_, "rhs vector of size " << rhs.size()
                                                       << " instead of " << n_);

    // initial guess
    Array_t<T> result = rhs;

    // solve tridiagonal system with SOR technique
    T omega = 1.5;
    T err = 2.0 * tol;
    T temp;
    for (Size sorIteration = 0; err > tol; ++sorIteration) {
        QL_REQUIRE(sorIteration < 100000, "tolerance ("
                                              << tol << ") not reached in "
                                              << sorIteration << " iterations. "
                                              << "The error still is " << err);

        temp = omega * (rhs[0] - upperDiagonal_[0] * result[1] -
                        diagonal_[0] * result[0]) /
               diagonal_[0];
        err = temp * temp;
        result[0] += temp;
        Size i;
        for (i = 1; i < n_ - 1; ++i) {
            temp = omega * (rhs[i] - upperDiagonal_[i] * result[i + 1] -
                            diagonal_[i] * result[i] -
                            lowerDiagonal_[i - 1] * result[i - 1]) /
                   diagonal_[i];
            err += temp * temp;
            result[i] += temp;
        }

        temp = omega * (rhs[i] - diagonal_[i] * result[i] -
                        lowerDiagonal_[i - 1] * result[i - 1]) /
               diagonal_[i];
        err += temp * temp;
        result[i] += temp;
    }
    return result;
}

template <class T>
Disposable<TridiagonalOperator_t<T> >
TridiagonalOperator_t<T>::identity(Size size) {
    TridiagonalOperator_t<T> I(Array_t<T>(size - 1, 0.0),  // lower diagonal
                               Array_t<T>(size, 1.0),      // diagonal
                               Array_t<T>(size - 1, 0.0)); // upper diagonal
    return I;
}

} // namespace QuantLib

#endif
