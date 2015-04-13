/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

/*! \file projectedconstraint.hpp
    \brief Projected constraint
*/

#ifndef quantlib_optimization_projectedconstraint_h
#define quantlib_optimization_projectedconstraint_h

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/projection.hpp>

namespace QuantLib {

template <class T> class ProjectedConstraint_t : public Constraint_t<T> {

  private:
    class Impl : public Constraint_t<T>::Impl {
      public:
        Impl(const Constraint_t<T> &constraint,
             const Array_t<T> &parameterValues,
             const std::vector<bool> &fixParameters)
            : constraint_(constraint),
              projection_(parameterValues, fixParameters) {}
        Impl(const Constraint_t<T> &constraint,
             const Projection_t<T> &projection)
            : constraint_(constraint), projection_(projection) {}
        bool test(const Array_t<T> &params) const {
            return constraint_.test(projection_.include(params));
        }
        Array_t<T> upperBound(const Array_t<T> &params) const {
            return constraint_.upperBound(projection_.include(params));
        }
        Array_t<T> lowerBound(const Array_t<T> &params) const {
            return constraint_.lowerBound(projection_.include(params));
        }

      private:
        const Constraint_t<T> &constraint_;
        const Projection_t<T> projection_;
    };

  public:
    ProjectedConstraint_t(const Constraint_t<T> &constraint,
                          const Array_t<T> &parameterValues,
                          const std::vector<bool> &fixParameters)
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new ProjectedConstraint_t<T>::Impl(constraint, parameterValues,
                                                 fixParameters))) {}

    ProjectedConstraint_t(const Constraint_t<T> &constraint,
                          const Projection_t<T> &projection)
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new ProjectedConstraint_t<T>::Impl(constraint, projection))) {}
};

typedef ProjectedConstraint_t<Real> ProjectedConstraint;

} // namespace QuantLib

#endif
