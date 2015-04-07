/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2005, 2007 StatPro Italia srl
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

/*! \file model.hpp
    \brief Abstract interest rate model class
*/

#ifndef quantlib_interest_rate_model_hpp
#define quantlib_interest_rate_model_hpp

#include <ql/option.hpp>
#include <ql/methods/lattices/lattice.hpp>
#include <ql/models/parameter.hpp>
#include <ql/models/calibrationhelper.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/projection.hpp>
#include <ql/math/optimization/projectedconstraint.hpp>

namespace QuantLib {

//! Affine model class
/*! Base class for analytically tractable models.

    \ingroup shortrate
*/
template <class T> class AffineModel_t : public virtual Observable {
  public:
    //! Implied discount curve
    virtual DiscountFactor discount(Time t) const = 0;

    virtual T discountBond(Time now, Time maturity,
                           Array_t<T> factors) const = 0;

    virtual T discountBondOption(Option::Type type, T strike, Time maturity,
                                 Time bondMaturity) const = 0;

    virtual T discountBondOption(Option::Type type, T strike, Time maturity,
                                 Time bondStart, Time bondMaturity) const {
        return discountBondOption(type, strike, maturity, bondMaturity);
    }
};

typedef AffineModel_t<Real> AffineModel;

//! Term-structure consistent model class
/*! This is a base class for models that can reprice exactly
    any discount bond.

    \ingroup shortrate
*/
template <class T>
class TermStructureConsistentModel_t : public virtual Observable {
  public:
    TermStructureConsistentModel_t(
        const Handle<YieldTermStructure_t<T> > &termStructure)
        : termStructure_(termStructure) {}
    const Handle<YieldTermStructure_t<T> > &termStructure() const {
        return termStructure_;
    }

  private:
    Handle<YieldTermStructure_t<T> > termStructure_;
};

typedef TermStructureConsistentModel_t<Real> TermStructureConsistentModel;

//! Calibrated model class
template <class T>
class CalibratedModel_t : public virtual Observer, public virtual Observable {
  public:
    CalibratedModel_t(Size nArguments);

    void update() {
        generateArguments();
        notifyObservers();
    }

    //! Calibrate to a set of market instruments (caps/swaptions)
    /*! An additional constraint can be passed which must be
        satisfied in addition to the constraints of the model.
    */
    virtual void
    calibrate(const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > > &,
              OptimizationMethod_t<T> &method,
              const EndCriteria_t<T> &endCriteria,
              const Constraint_t<T> &constraint = Constraint_t<T>(),
              const std::vector<T> &weights = std::vector<T>(),
              const std::vector<bool> &fixParameters = std::vector<bool>());

    T value(const Array_t<T> &params,
            const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > > &);

    const boost::shared_ptr<Constraint_t<T> > &constraint() const;
    //! returns end criteria result
    typename EndCriteria_t<T>::Type endCriteria();
    //! Returns array of arguments on which calibration is done
    Disposable<Array_t<T> > params() const;

    virtual void setParams(const Array_t<T> &params);

  protected:
    virtual void generateArguments() {}
    std::vector<Parameter> arguments_;
    boost::shared_ptr<Constraint_t<T> > constraint_;
    typename EndCriteria_t<T>::Type shortRateEndCriteria_;

  private:
    //! Constraint_t<T> imposed on arguments
    class PrivateConstraint;
    //! Calibration cost function class
    class CalibrationFunction;
    friend class CalibrationFunction;
};

typedef CalibratedModel_t<Real> CalibratedModel;

//! Abstract short-rate model class
/*! \ingroup shortrate */
template <class T> class ShortRateModel_t : public CalibratedModel_t<T> {
  public:
    ShortRateModel_t(Size nArguments);
    virtual boost::shared_ptr<Lattice> tree(const TimeGrid &) const = 0;
};

typedef ShortRateModel_t<Real> ShortRateModel;

// inline definitions

template <class T>
inline const boost::shared_ptr<Constraint_t<T> > &
CalibratedModel_t<T>::constraint() const {
    return constraint_;
}

template <class T>
class CalibratedModel_t<T>::PrivateConstraint : public Constraint_t<T> {
  private:
    class Impl : public Constraint_t<T>::Impl {
      public:
        Impl(const std::vector<Parameter> &arguments) : arguments_(arguments) {}

        bool test(const Array_t<T> &params) const {
            Size k = 0;
            for (Size i = 0; i < arguments_.size(); i++) {
                Size size = arguments_[i].size();
                Array_t<T> testParams(size);
                for (Size j = 0; j < size; j++, k++)
                    testParams[j] = params[k];
                if (!arguments_[i].testParams(testParams))
                    return false;
            }
            return true;
        }

        Array_t<T> upperBound(const Array_t<T> &params) const {
            Size k = 0, k2 = 0;
            Size totalSize = 0;
            for (Size i = 0; i < arguments_.size(); i++) {
                totalSize += arguments_[i].size();
            }
            Array_t<T> result(totalSize);
            for (Size i = 0; i < arguments_.size(); i++) {
                Size size = arguments_[i].size();
                Array_t<T> partialParams(size);
                for (Size j = 0; j < size; j++, k++)
                    partialParams[j] = params[k];
                Array_t<T> tmpBound =
                    arguments_[i].constraint().upperBound(partialParams);
                for (Size j = 0; j < size; j++, k2++)
                    result[k2] = tmpBound[j];
            }
            return result;
        }

        Array_t<T> lowerBound(const Array_t<T> &params) const {
            Size k = 0, k2 = 0;
            Size totalSize = 0;
            for (Size i = 0; i < arguments_.size(); i++) {
                totalSize += arguments_[i].size();
            }
            Array_t<T> result(totalSize);
            for (Size i = 0; i < arguments_.size(); i++) {
                Size size = arguments_[i].size();
                Array_t<T> partialParams(size);
                for (Size j = 0; j < size; j++, k++)
                    partialParams[j] = params[k];
                Array_t<T> tmpBound =
                    arguments_[i].constraint().lowerBound(partialParams);
                for (Size j = 0; j < size; j++, k2++)
                    result[k2] = tmpBound[j];
            }
            return result;
        }

      private:
        const std::vector<Parameter> &arguments_;
    };

  public:
    PrivateConstraint(const std::vector<Parameter_t<T> > &arguments)
        : Constraint_t<T>(boost::shared_ptr<typename Constraint_t<T>::Impl>(
              new PrivateConstraint::Impl(arguments))) {}
};

// implementation

template <class T>
CalibratedModel_t<T>::CalibratedModel_t(Size nArguments)
    : arguments_(nArguments), constraint_(new PrivateConstraint(arguments_)),
      shortRateEndCriteria_(EndCriteria_t<T>::None) {}

template <class T>
class CalibratedModel_t<T>::CalibrationFunction : public CostFunction_t<T> {
  public:
    CalibrationFunction(
        CalibratedModel_t<T> *model,
        const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > >
            &instruments,
        const std::vector<T> &weights, const Projection_t<T> &projection)
        : model_(model, no_deletion), instruments_(instruments),
          weights_(weights), projection_(projection) {}

    virtual ~CalibrationFunction() {}

    virtual T value(const Array_t<T> &params) const {
        model_->setParams(projection_.include(params));
        T value = 0.0;
        for (Size i = 0; i < instruments_.size(); i++) {
            T diff = instruments_[i]->calibrationError();
            value += diff * diff * weights_[i];
        }
        return QLFCT::sqrt(value);
    }

    virtual Disposable<Array_t<T> > values(const Array_t<T> &params) const {
        model_->setParams(projection_.include(params));
        Array_t<T> values(instruments_.size());
        for (Size i = 0; i < instruments_.size(); i++) {
            values[i] =
                instruments_[i]->calibrationError() * QLFCT::sqrt(weights_[i]);
        }
        return values;
    }

    virtual T finiteDifferenceEpsilon() const { return 1e-6; }

  private:
    boost::shared_ptr<CalibratedModel_t<T> > model_;
    const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > >
        &instruments_;
    std::vector<T> weights_;
    const Projection_t<T> projection_;
};

template <class T>
void CalibratedModel_t<T>::calibrate(
    const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > > &instruments,
    OptimizationMethod_t<T> &method, const EndCriteria_t<T> &endCriteria,
    const Constraint_t<T> &additionalConstraint, const std::vector<T> &weights,
    const std::vector<bool> &fixParameters) {

    QL_REQUIRE(weights.empty() || weights.size() == instruments.size(),
               "mismatch between number of instruments and weights");

    Constraint_t<T> c;
    if (additionalConstraint.empty())
        c = *constraint_;
    else
        c = CompositeConstraint(*constraint_, additionalConstraint);
    std::vector<T> w =
        weights.empty() ? std::vector<T>(instruments.size(), 1.0) : weights;

    Array_t<T> prms = params();
    std::vector<bool> all(prms.size(), false);
    Projection_t<T> proj(prms, fixParameters.size() > 0 ? fixParameters : all);
    CalibrationFunction f(this, instruments, w, proj);
    ProjectedConstraint_t<T> pc(c, proj);
    Problem_t<T> prob(f, pc, proj.project(prms));
    shortRateEndCriteria_ = method.minimize(prob, endCriteria);
    Array_t<T> result(prob.currentValue());
    setParams(proj.include(result));
    Array_t<T> shortRateProblemValues_ = prob.values(result);

    notifyObservers();
}

template <class T>
typename EndCriteria_t<T>::Type CalibratedModel_t<T>::endCriteria() {
    return shortRateEndCriteria_;
}

template <class T>
T CalibratedModel_t<T>::value(
    const Array_t<T> &params,
    const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > >
        &instruments) {
    std::vector<T> w = std::vector<T>(instruments.size(), 1.0);
    Projection_t<T> p(params);
    CalibrationFunction f(this, instruments, w, p);
    return f.value(params);
}

template <class T>
Disposable<Array_t<T> > CalibratedModel_t<T>::params() const {
    Size size = 0, i;
    for (i = 0; i < arguments_.size(); i++)
        size += arguments_[i].size();
    Array_t<T> params(size);
    Size k = 0;
    for (i = 0; i < arguments_.size(); i++) {
        for (Size j = 0; j < arguments_[i].size(); j++, k++) {
            params[k] = arguments_[i].params()[j];
        }
    }
    return params;
}

template <class T>
void CalibratedModel_t<T>::setParams(const Array_t<T> &params) {
    typename Array_t<T>::const_iterator p = params.begin();
    for (Size i = 0; i < arguments_.size(); ++i) {
        for (Size j = 0; j < arguments_[i].size(); ++j, ++p) {
            QL_REQUIRE(p != params.end(), "parameter array too small");
            arguments_[i].setParam(j, *p);
        }
    }
    QL_REQUIRE(p == params.end(), "parameter array too big!");
    generateArguments();
    notifyObservers();
}

template <class T>
ShortRateModel_t<T>::ShortRateModel_t(Size nArguments)
    : CalibratedModel_t<T>(nArguments) {}
} // namespace quantlib

#endif
