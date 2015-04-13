/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2012 Mateusz Kapturski
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

/*! \file constraint.hpp
    \brief Abstract constraint class
*/

#ifndef quantlib_optimization_constraint_h
#define quantlib_optimization_constraint_h

#include <ql/math/array.hpp>

namespace QuantLib {

//! Base constraint class
template <class T> class Constraint_t {
  protected:
    //! Base class for constraint implementations
    class Impl {
      public:
        virtual ~Impl() {}
        //! Tests if params satisfy the constraint
        virtual bool test(const Array_t<T> &params) const = 0;
        //! Returns upper bound for given parameters
        virtual Array_t<T> upperBound(const Array_t<T> &params) const {
            return Array_t<T>(params.size(),
                              std::numeric_limits<Array::value_type>::max());
        }
        //! Returns lower bound for given parameters
        virtual Array_t<T> lowerBound(const Array_t<T> &params) const {
            return Array_t<T>(params.size(),
                              std::numeric_limits<Array::value_type>::min());
        }
    };
    boost::shared_ptr<Impl> impl_;

  public:
    bool empty() const { return !impl_; }
    bool test(const Array_t<T> &p) const { return impl_->test(p); }
    Array_t<T> upperBound(const Array_t<T> &params) const {
        Array_t<T> result = impl_->upperBound(params);
        QL_REQUIRE(params.size() == result.size(),
                   "upper bound size (" << result.size()
                                        << ") not equal to params size ("
                                        << params.size() << ")");
        return result;
    }
    Array_t<T> lowerBound(const Array_t<T> &params) const {
        Array_t<T> result = impl_->lowerBound(params);
        QL_REQUIRE(params.size() == result.size(),
                   "lower bound size (" << result.size()
                                        << ") not equal to params size ("
                                        << params.size() << ")");
        return result;
    }
    T update(Array_t<T> &p, const Array_t<T> &direction, T beta);
    Constraint_t(
        const boost::shared_ptr<Impl> &impl = boost::shared_ptr<Impl>());
};

typedef Constraint_t<Real> Constraint;

//! No constraint
template <class T> class NoConstraint_t : public Constraint_t<T> {
  private:
    class Impl : public Constraint_t<T>::Impl {
      public:
        bool test(const Array_t<T> &) const { return true; }
    };

  public:
    NoConstraint_t()
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new NoConstraint_t<T>::Impl)) {}
};

typedef NoConstraint_t<Real> NoConstraint;

//! %Constraint_t imposing positivity to all arguments
template <class T> class PositiveConstraint_t : public Constraint_t<T> {
  private:
    class Impl : public Constraint_t<T>::Impl {
      public:
        bool test(const Array_t<T> &params) const {
            for (Size i = 0; i < params.size(); ++i) {
                if (params[i] <= 0.0)
                    return false;
            }
            return true;
        }
        Array_t<T> upperBound(const Array_t<T> &params) const {
            return Array_t<T>(params.size(),
                              std::numeric_limits<Array::value_type>::max());
        }
        Array_t<T> lowerBound(const Array_t<T> &params) const {
            return Array_t<T>(params.size(), 0.0);
        }
    };

  public:
    PositiveConstraint_t()
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new PositiveConstraint_t<T>::Impl)) {}
};

typedef PositiveConstraint_t<Real> PositiveConstraint;

//! %Constraint_t imposing all arguments to be in [low,high]
template <class T> class BoundaryConstraint_t : public Constraint_t<T> {
  private:
    class Impl : public Constraint_t<T>::Impl {
      public:
        Impl(T low, T high) : low_(low), high_(high) {}
        bool test(const Array_t<T> &params) const {
            for (Size i = 0; i < params.size(); i++) {
                if ((params[i] < low_) || (params[i] > high_))
                    return false;
            }
            return true;
        }
        Array_t<T> upperBound(const Array_t<T> &params) const {
            return Array_t<T>(params.size(), high_);
        }
        Array_t<T> lowerBound(const Array_t<T> &params) const {
            return Array_t<T>(params.size(), low_);
        }

      private:
        T low_, high_;
    };

  public:
    BoundaryConstraint_t(T low, T high)
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new BoundaryConstraint_t<T>::Impl(low, high))) {}
};

typedef BoundaryConstraint_t<Real> BoundaryConstraint;

//! %Constraint enforcing both given sub-constraints
template <class T> class CompositeConstraint_t : public Constraint_t<T> {
  private:
    class Impl : public Constraint_t<T>::Impl {
      public:
        Impl(const Constraint_t<T> &c1, const Constraint_t<T> &c2)
            : c1_(c1), c2_(c2) {}
        bool test(const Array_t<T> &params) const {
            return c1_.test(params) && c2_.test(params);
        }
        Array_t<T> upperBound(const Array_t<T> &params) const {
            Array_t<T> c1ub = c1_.upperBound(params);
            Array_t<T> c2ub = c2_.upperBound(params);
            Array_t<T> rtrnArray(c1ub.size(), 0.0);
            for (Size iter = 0; iter < c1ub.size(); iter++) {
                rtrnArray.at(iter) = QLFCT::min(c1ub.at(iter), c2ub.at(iter));
            }
            return rtrnArray;
        }
        Array_t<T> lowerBound(const Array_t<T> &params) const {
            Array_t<T> c1lb = c1_.lowerBound(params);
            Array_t<T> c2lb = c2_.lowerBound(params);
            Array_t<T> rtrnArray(c1lb.size(), 0.0);
            for (Size iter = 0; iter < c1lb.size(); iter++) {
                rtrnArray.at(iter) = QLFCT::max(c1lb.at(iter), c2lb.at(iter));
            }
            return rtrnArray;
        }

      private:
        Constraint_t<T> c1_, c2_;
    };

  public:
    CompositeConstraint_t(const Constraint_t<T> &c1, const Constraint_t<T> &c2)
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new CompositeConstraint_t<T>::Impl(c1, c2))) {}
};

typedef CompositeConstraint_t<Real> CompositeConstraint;

//! %Constraint imposing i-th argument to be in [low_i,high_i] for all i
template <class T>
class NonhomogeneousBoundaryConstraint_t : public Constraint_t<T> {
  private:
    class Impl : public Constraint_t<T>::Impl {
      public:
        Impl(const Array_t<T> &low, const Array_t<T> &high)
            : low_(low), high_(high) {
            QL_ENSURE(low_.size() == high_.size(),
                      "Upper and lower boundaries sizes are inconsistent.");
        }
        bool test(const Array_t<T> &params) const {
            QL_ENSURE(
                params.size() == low_.size(),
                "Number of parameters and boundaries sizes are inconsistent.")
            for (Size i = 0; i < params.size(); i++) {
                if ((params[i] < low_[i]) || (params[i] > high_[i]))
                    return false;
            }
            return true;
        }
        Array_t<T> upperBound(const Array_t<T> &) const { return high_; }
        Array_t<T> lowerBound(const Array_t<T> &) const { return low_; }

      private:
        Array_t<T> low_, high_;
    };

  public:
    NonhomogeneousBoundaryConstraint_t(Array low, Array_t<T> high)
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new NonhomogeneousBoundaryConstraint_t<T>::Impl(low, high))) {}
};

typedef NonhomogeneousBoundaryConstraint_t<Real>
    NonhomogeneousBoundaryConstraint;

// implementation

template <class T>
Constraint_t<T>::Constraint_t(
    const boost::shared_ptr<Constraint_t<T>::Impl> &impl)
    : impl_(impl) {}

template <class T>
T Constraint_t<T>::update(Array_t<T> &params, const Array_t<T> &direction,
                          T beta) {

    T diff = beta;
    Array_t<T> newParams = params + diff * direction;
    bool valid = test(newParams);
    Integer icount = 0;
    while (!valid) {
        if (icount > 200)
            QL_FAIL("can't update parameter vector");
        diff *= 0.5;
        icount++;
        newParams = params + diff * direction;
        valid = test(newParams);
    }
    params += diff * direction;
    return diff;
}
}

#endif
