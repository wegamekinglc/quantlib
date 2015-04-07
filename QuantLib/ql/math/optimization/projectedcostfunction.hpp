/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 François du Vignaud
 Copyright (C) 2007 Giorgio Facchinetti
 Copyright (C) 2013, 2015 Peter Caspers

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

/*! \file projectedcostfunction.hpp
    \brief Cost function utility
*/

#ifndef quantlib_math_projectedcostfunction_h
#define quantlib_math_projectedcostfunction_h

#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/projection.hpp>

namespace QuantLib {

//! Parameterized cost function
/*! This class creates a proxy cost function which can depend
    on any arbitrary subset of parameters (the other being fixed)
*/

template <class T>
class ProjectedCostFunction_t : public CostFunction_t<T>,
                                public Projection_t<T> {
  public:
    ProjectedCostFunction_t(const CostFunction &costFunction,
                            const Array_t<T> &parameterValues,
                            const std::vector<bool> &fixParameters);

    ProjectedCostFunction_t(const CostFunction &costFunction,
                            const Projection &projection);

    //! \name CostFunction interface
    //@{
    virtual T value(const Array_t<T> &freeParameters) const;
    virtual Disposable<Array_t<T> >
    values(const Array_t<T> &freeParameters) const;
    //@}

  private:
    const CostFunction_t<T> &costFunction_;
};

typedef ProjectedCostFunction_t<Real> ProjectedCostFunction;

// implementation

template <class T>
ProjectedCostFunction_t<T>::ProjectedCostFunction_t(
    const CostFunction &costFunction, const Array_t<T> &parameterValues,
    const std::vector<bool> &fixParameters)
    : Projection(parameterValues, fixParameters), costFunction_(costFunction) {}

template <class T>
ProjectedCostFunction_t<T>::ProjectedCostFunction_t(
    const CostFunction &costFunction, const Projection &projection)
    : Projection(projection), costFunction_(costFunction) {}

template <class T>
T ProjectedCostFunction_t<T>::value(const Array_t<T> &freeParameters) const {
    this->mapFreeParameters(freeParameters);
    return costFunction_.value(this->actualParameters_);
}

template <class T>
Disposable<Array_t<T> >
ProjectedCostFunction_t<T>::values(const Array_t<T> &freeParameters) const {
    this->mapFreeParameters(freeParameters);
    return costFunction_.values(this->actualParameters_);
}
}

#endif
