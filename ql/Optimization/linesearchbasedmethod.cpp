/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/reference/license.html>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/Optimization/linesearchbasedmethod.hpp>
#include <ql/Optimization/armijo.hpp>

namespace QuantLib {

    LineSearchBasedMethod::LineSearchBasedMethod(
                           const boost::shared_ptr<LineSearch>& lineSearch)
    : lineSearch_(lineSearch) {
        if (!lineSearch_)
           lineSearch_ = boost::shared_ptr<LineSearch>(new ArmijoLineSearch);
    }

}
