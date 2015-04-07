/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
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

/*! \file parameter.hpp
    \brief Model parameter classes
*/

#ifndef quantlib_interest_rate_modelling_parameter_hpp
#define quantlib_interest_rate_modelling_parameter_hpp

#include <ql/qldefines.hpp>
#include <ql/handle.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <vector>

namespace QuantLib {

//! Base class for model arguments
template <class T> class Parameter_t {
  protected:
    //! Base class for model parameter implementation
    class Impl {
      public:
        virtual ~Impl() {}
        virtual T value(const Array_t<T> &params, Time t) const = 0;
    };
    boost::shared_ptr<Impl> impl_;

  public:
    Parameter_t() : constraint_(NoConstraint()) {}
    const Array_t<T> &params() const { return params_; }
    void setParam(Size i, T x) { params_[i] = x; }
    bool testParams(const Array_t<T> &params) const {
        return constraint_.test(params);
    }
    Size size() const { return params_.size(); }
    T operator()(Time t) const { return impl_->value(params_, t); }
    const boost::shared_ptr<Impl> &implementation() const { return impl_; }
    const Constraint_t<T> &constraint() const { return constraint_; }

  protected:
    Parameter_t(Size size, const boost::shared_ptr<Impl> &impl,
                const Constraint_t<T> &constraint)
        : impl_(impl), params_(size), constraint_(constraint) {}
    Array_t<T> params_;
    Constraint_t<T> constraint_;
};

typedef Parameter_t<Real> Parameter;

//! Standard constant parameter \f$ a(t) = a \f$
template <class T> class ConstantParameter_t : public Parameter_t<T> {
  private:
    class Impl : public Parameter_t<T>::Impl {
      public:
        T value(const Array_t<T> &params, Time) const { return params[0]; }
    };

  public:
    ConstantParameter_t(const Constraint_t<T> &constraint)
        : Parameter_t<T>(1, boost::shared_ptr<typename Parameter_t<T>::Impl>(
                             new ConstantParameter_t<T>::Impl),
                      constraint) {}

    ConstantParameter_t(T value, const Constraint_t<T> &constraint)
        : Parameter_t<T>(1, boost::shared_ptr<typename Parameter_t<T>::Impl>(
                             new ConstantParameter_t<T>::Impl),
                      constraint) {
        this->params_[0] = value;
        QL_REQUIRE(this->testParams(this->params_), value << ": invalid value");
    }
};

typedef ConstantParameter_t<Real> ConstantParameter;

//! %Parameter which is always zero \f$ a(t) = 0 \f$
template <class T> class NullParameter_t : public Parameter_t<T> {
  private:
    class Impl : public Parameter_t<T>::Impl {
      public:
        T value(const Array_t<T> &, Time) const { return 0.0; }
    };

  public:
    NullParameter_t()
        : Parameter_t<T>(0, boost::shared_ptr<typename Parameter_t<T>::Impl>(
                             new NullParameter_t<T>::Impl),
                      NoConstraint()) {}
};

typedef NullParameter_t<Real> NullParameter;

//! Piecewise-constant parameter
/*! \f$ a(t) = a_i if t_{i-1} \geq t < t_i \f$.
    This kind of parameter is usually used to enhance the fitting of a
    model
*/
template <class T> class PiecewiseConstantParameter_t : public Parameter_t<T> {
  private:
    class Impl : public Parameter_t<T>::Impl {
      public:
        Impl(const std::vector<Time> &times) : times_(times) {}

        T value(const Array_t<T> &params, Time t) const {
            Size size = times_.size();
            for (Size i = 0; i < size; i++) {
                if (t < times_[i])
                    return params[i];
            }
            return params[size];
        }

      private:
        std::vector<Time> times_;
    };

  public:
    PiecewiseConstantParameter_t(
        const std::vector<Time> &times,
        const Constraint_t<T> &constraint = NoConstraint())
        : Parameter_t<T>(times.size() + 1,
                      boost::shared_ptr<typename Parameter_t<T>::Impl>(
                          new PiecewiseConstantParameter_t<T>::Impl(times)),
                      constraint) {}
};

typedef PiecewiseConstantParameter_t<Real> PiecewiseConstantParameter;

//! Deterministic time-dependent parameter used for yield-curve fitting
template <class T>
class TermStructureFittingParameter_t : public Parameter_t<T> {
  public:
    class NumericalImpl : public Parameter_t<T>::Impl {
      public:
        NumericalImpl(const Handle<YieldTermStructure_t<T> > &termStructure)
            : times_(0), values_(0), termStructure_(termStructure) {}

        void set(Time t, T x) {
            times_.push_back(t);
            values_.push_back(x);
        }
        void change(T x) { values_.back() = x; }
        void reset() {
            times_.clear();
            values_.clear();
        }
        T value(const Array_t<T> &, Time t) const {
            std::vector<Time>::const_iterator result =
                std::find(times_.begin(), times_.end(), t);
            QL_REQUIRE(result != times_.end(), "fitting parameter not set!");
            return values_[result - times_.begin()];
        }
        const Handle<YieldTermStructure_t<T> > &termStructure() const {
            return termStructure_;
        }

      private:
        std::vector<Time> times_;
        std::vector<T> values_;
        Handle<YieldTermStructure_t<T> > termStructure_;
    };

    TermStructureFittingParameter_t(
        const boost::shared_ptr<typename Parameter_t<T>::Impl> &impl)
        : Parameter_t<T>(0, impl, NoConstraint()) {}

    TermStructureFittingParameter_t(
        const Handle<YieldTermStructure_t<T> > &term)
        : Parameter_t<T>(0, boost::shared_ptr<typename Parameter_t<T>::Impl>(
                             new NumericalImpl(term)),
                      NoConstraint()) {}
};

typedef TermStructureFittingParameter_t<Real> TermStructureFittingParameter;

} // namespace QuantLib

#endif
