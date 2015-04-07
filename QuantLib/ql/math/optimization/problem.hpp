/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 Ferdinando Ametrano
 Copyright (C) 2007 François du Vignaud
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

/*! \file problem.hpp
    \brief Abstract optimization problem class
*/

#ifndef quantlib_optimization_problem_h
#define quantlib_optimization_problem_h

#include <ql/math/optimization/method.hpp>
#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/constraint.hpp>

namespace QuantLib {

//! Constrained optimization problem
/*! \warning The passed CostFunction and Constraint instances are
             stored by reference.  The user of this class must
             make sure that they are not destroyed before the
             Problem instance.
*/
template <class T> class Problem_t {
  public:
    //! default constructor
    Problem_t(CostFunction_t<T> &costFunction, Constraint_t<T> &constraint,
              const Array_t<T> &initialValue = Array_t<T>())
        : costFunction_(costFunction), constraint_(constraint),
          currentValue_(initialValue) {}

    /*! \warning it does not reset the current minumum to any initial value
    */
    void reset();

    //! call cost function computation and increment evaluation counter
    T value(const Array_t<T> &x);

    //! call cost values computation and increment evaluation counter
    Disposable<Array_t<T> > values(const Array_t<T> &x);

    //! call cost function gradient computation and increment
    //  evaluation counter
    void gradient(Array_t<T> &grad_f, const Array_t<T> &x);

    //! call cost function computation and it gradient
    T valueAndGradient(Array_t<T> &grad_f, const Array_t<T> &x);

    //! Constraint
    Constraint_t<T> &constraint() const { return constraint_; }

    //! Cost function
    CostFunction_t<T> &costFunction() const { return costFunction_; }

    void setCurrentValue(const Array_t<T> &currentValue) {
        currentValue_ = currentValue;
    }

    //! current value of the local minimum
    const Array_t<T> &currentValue() const { return currentValue_; }

    void setFunctionValue(T functionValue) { functionValue_ = functionValue; }

    //! value of cost function
    T functionValue() const { return functionValue_; }

    void setGradientNormValue(T squaredNorm) { squaredNorm_ = squaredNorm; }
    //! value of cost function gradient norm
    T gradientNormValue() const { return squaredNorm_; }

    //! number of evaluation of cost function
    Integer functionEvaluation() const { return functionEvaluation_; }

    //! number of evaluation of cost function gradient
    Integer gradientEvaluation() const { return gradientEvaluation_; }

  protected:
    //! Unconstrained cost function
    CostFunction_t<T> &costFunction_;
    //! Constraint
    Constraint_t<T> &constraint_;
    //! current value of the local minimum
    Array_t<T> currentValue_;
    //! function and gradient norm values at the curentValue_ (i.e. the last
    // step)
    T functionValue_, squaredNorm_;
    //! number of evaluation of cost function and its gradient
    Integer functionEvaluation_, gradientEvaluation_;
};

typedef Problem_t<Real> Problem;

// inline definitions
template <class T> inline T Problem_t<T>::value(const Array_t<T> &x) {
    ++functionEvaluation_;
    return costFunction_.value(x);
}

template <class T>
inline Disposable<Array_t<T> > Problem_t<T>::values(const Array_t<T> &x) {
    ++functionEvaluation_;
    return costFunction_.values(x);
}

template <class T>
inline void Problem_t<T>::gradient(Array_t<T> &grad_f, const Array_t<T> &x) {
    ++gradientEvaluation_;
    costFunction_.gradient(grad_f, x);
}

template <class T>
inline T Problem_t<T>::valueAndGradient(Array_t<T> &grad_f,
                                        const Array_t<T> &x) {
    ++functionEvaluation_;
    ++gradientEvaluation_;
    return costFunction_.valueAndGradient(grad_f, x);
}

template <class T> inline void Problem_t<T>::reset() {
    functionEvaluation_ = gradientEvaluation_ = 0;
    functionValue_ = squaredNorm_ = Null<T>();
}
}

#endif
