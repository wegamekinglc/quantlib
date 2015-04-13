/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2001, 2002, 2003 Nicolas Di Césaré
 Copyright (C) 2007 François du Vignaud
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

/*! \file method.hpp
    \brief Abstract optimization method class
*/

#ifndef quantlib_optimization_method_h
#define quantlib_optimization_method_h

#include <ql/math/optimization/endcriteria.hpp>

namespace QuantLib {

template <class T> class Problem_t;

//! Abstract class for constrained optimization method
template <class T> class OptimizationMethod_t {
  public:
    virtual ~OptimizationMethod_t() {}

    //! minimize the optimization problem P
    virtual typename EndCriteria_t<T>::Type
    minimize(Problem_t<T> &P, const EndCriteria_t<T> &endCriteria) = 0;
};

typedef OptimizationMethod_t<Real> OptimizationMethod;
}

#endif
