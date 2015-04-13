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

/*! \file gsr.hpp
    \brief GSR 1 factor model
*/

#ifndef quantlib_gsr_hpp
#define quantlib_gsr_hpp

#include <ql/time/schedule.hpp>
#include <ql/math/integrals/simpsonintegral.hpp>
#include <ql/math/integrals/gausslobattointegral.hpp>
#include <ql/math/distributions/normaldistribution.hpp>

#include <ql/experimental/models/gaussian1dmodel.hpp>
#include <ql/experimental/models/gsrprocess.hpp>
#include <ql/quotes/simplequote.hpp>

#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>

namespace QuantLib {

//! One factor gsr model, formulation is in forward measure

template <class T>
class Gsr_t : public Gaussian1dModel_t<T>, public CalibratedModel_t<T> {

  public:
    // constant mean reversion
    Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
          const std::vector<Date> &volstepdates,
          const std::vector<T> &volatilities, const T reversion,
          const Time T0 = 60.0);
    // piecewise mean reversion (with same step dates as volatilities)
    Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
          const std::vector<Date> &volstepdates,
          const std::vector<T> &volatilities, const std::vector<T> &reversions,
          const Time T0 = 60.0);
    // constant mean reversion with floating model data
    Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
          const std::vector<Date> &volstepdates,
          const std::vector<Handle<Quote_t<T> > > &volatilities,
          const Handle<Quote_t<T> > reversion, const Time T0 = 60.0);
    // piecewise mean reversion with floating model data
    Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
          const std::vector<Date> &volstepdates,
          const std::vector<Handle<Quote_t<T> > > &volatilities,
          const std::vector<Handle<Quote_t<T> > > &reversions,
          const Time T0 = 60.0);

    const T numeraireTime() const;
    const void numeraireTime(const Time T0);

    const Array_t<T> &reversion() const { return reversion_.params(); }
    const Array_t<T> &volatility() const { return sigma_.params(); }

    // calibration constraints

    Disposable<std::vector<bool> > FixedReversions() {
        std::vector<bool> res(reversions_.size(), true);
        std::vector<bool> vol(volatilities_.size(), false);
        res.insert(res.end(), vol.begin(), vol.end());
        return res;
    }

    Disposable<std::vector<bool> > FixedVolatilities() {
        std::vector<bool> res(reversions_.size(), false);
        std::vector<bool> vol(volatilities_.size(), true);
        res.insert(res.end(), vol.begin(), vol.end());
        return res;
    }

    Disposable<std::vector<bool> > MoveVolatility(Size i) {
        QL_REQUIRE(i < volatilities_.size(),
                   "volatility with index " << i << " does not exist (0..."
                                            << volatilities_.size() - 1 << ")");
        std::vector<bool> res(reversions_.size() + volatilities_.size(), true);
        res[reversions_.size() + i] = false;
        return res;
    }

    Disposable<std::vector<bool> > MoveReversion(Size i) {
        QL_REQUIRE(i < reversions_.size(),
                   "reversion with index " << i << " does not exist (0..."
                                           << reversions_.size() - 1 << ")");
        std::vector<bool> res(reversions_.size() + volatilities_.size(), true);
        res[i] = false;
        return res;
    }

    // With fixed reversion calibrate the volatilities one by one
    // to the given helpers. It is assumed that that volatility step
    // dates are suitable for this, i.e. they should be identical to
    // the fixing dates of the helpers (except for the last one where
    // we do not need a step). Also note that the endcritera reflect
    // only the status of the last calibration when using this method.
    void calibrateVolatilitiesIterative(
        const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > > &helpers,
        OptimizationMethod_t<T> &method, const EndCriteria_t<T> &endCriteria,
        const Constraint_t<T> &constraint = Constraint_t<T>(),
        const std::vector<T> &weights = std::vector<T>()) {

        for (Size i = 0; i < helpers.size(); i++) {
            std::vector<boost::shared_ptr<CalibrationHelper_t<T> > > h(
                1, helpers[i]);
            this->calibrate(h, method, endCriteria, constraint, weights,
                            MoveVolatility(i));
        }
    }

    // With fixed volatility calibrate the reversions one by one
    // to the given helpers. In this case the step dates must be chosen
    // according to the maturities of the calibration instruments.
    void calibrateReversionsIterative(
        const std::vector<boost::shared_ptr<CalibrationHelper_t<T> > > &helpers,
        OptimizationMethod_t<T> &method, const EndCriteria_t<T> &endCriteria,
        const Constraint_t<T> &constraint = Constraint_t<T>(),
        const std::vector<T> &weights = std::vector<T>()) {

        for (Size i = 0; i < helpers.size(); i++) {
            std::vector<boost::shared_ptr<CalibrationHelper_t<T> > > h(
                1, helpers[i]);
            calibrate(h, method, endCriteria, constraint, weights,
                      MoveReversion(i));
        }
    }

  protected:
    const T numeraireImpl(const Time t, const T y,
                          const Handle<YieldTermStructure_t<T> > &yts) const;

    const T zerobondImpl(const Time T0, const Time t, const T y,
                         const Handle<YieldTermStructure_t<T> > &yts) const;

    void generateArguments() {
        boost::static_pointer_cast<GsrProcess_t<T> >(this->stateProcess_)
            ->flushCache();
        this->notifyObservers();
    }

    void update() { LazyObject::update(); }

    void performCalculations() const {
        Gaussian1dModel_t<T>::performCalculations();
        updateTimes();
        updateState();
    }

  private:
    void updateTimes() const;
    void updateState() const;
    void initialize(Time);

    Parameter_t<T> &reversion_, &sigma_;

    std::vector<Handle<Quote_t<T> > > volatilities_;
    std::vector<Handle<Quote_t<T> > > reversions_;
    std::vector<Date> volstepdates_; // this is shared between vols and
                                     // reverisons in case of piecewise
                                     // reversions
    mutable std::vector<Time> volsteptimes_;
    mutable Array_t<Time>
        volsteptimesArray_; // FIXME this is redundant (just a copy of
                            // volsteptimes_)
};

typedef Gsr_t<Real> Gsr;

// inline

template <class T> inline const T Gsr_t<T>::numeraireTime() const {
    return boost::dynamic_pointer_cast<GsrProcess_t<T> >(this->stateProcess_)
        ->getForwardMeasureTime();
}

template <class T> inline const void Gsr_t<T>::numeraireTime(const Time T0) {
    boost::dynamic_pointer_cast<GsrProcess_t<T> >(this->stateProcess_)
        ->setForwardMeasureTime(T0);
}

// implementation

template <class T>
Gsr_t<T>::Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
                const std::vector<Date> &volstepdates,
                const std::vector<T> &volatilities, const T reversion,
                const Time T0)
    : Gaussian1dModel_t<T>(termStructure), CalibratedModel_t<T>(2),
      reversion_(this->arguments_[0]), sigma_(this->arguments_[1]),
      volstepdates_(volstepdates) {

    QL_REQUIRE(!termStructure.empty(), "yield term structure handle is empty");

    volatilities_.resize(volatilities.size());
    for (Size i = 0; i < volatilities.size(); ++i)
        volatilities_[i] = Handle<Quote_t<T> >(
            boost::make_shared<SimpleQuote_t<T> >(volatilities[i]));
    reversions_.resize(1);
    reversions_[0] =
        Handle<Quote_t<T> >(boost::make_shared<SimpleQuote_t<T> >(reversion));

    initialize(T0);
}

template <class T>
Gsr_t<T>::Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
                const std::vector<Date> &volstepdates,
                const std::vector<T> &volatilities,
                const std::vector<T> &reversions, const Time T0)
    : Gaussian1dModel(termStructure), CalibratedModel_t<T>(2),
      reversion_(this->arguments_[0]), sigma_(this->arguments_[1]),
      volstepdates_(volstepdates) {

    QL_REQUIRE(!termStructure.empty(), "yield term structure handle is empty");

    volatilities_.resize(volatilities.size());
    for (Size i = 0; i < volatilities.size(); ++i)
        volatilities_[i] = Handle<Quote_t<T> >(
            boost::make_shared<SimpleQuote_t<T> >(volatilities[i]));
    reversions_.resize(reversions.size());
    for (Size i = 0; i < reversions.size(); ++i)
        reversions_[i] = Handle<Quote_t<T> >(
            boost::make_shared<SimpleQuote_t<T> >(reversions[i]));

    initialize(T0);
}

template <class T>
Gsr_t<T>::Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
                const std::vector<Date> &volstepdates,
                const std::vector<Handle<Quote_t<T> > > &volatilities,
                const Handle<Quote_t<T> > reversion, const Time T0)
    : Gaussian1dModel(termStructure), CalibratedModel_t<T>(2),
      reversion_(this->arguments_[0]), sigma_(this->arguments_[1]),
      volatilities_(volatilities),
      reversions_(std::vector<Handle<Quote_t<T> > >(1, reversion)),
      volstepdates_(volstepdates) {

    QL_REQUIRE(!termStructure.empty(), "yield term structure handle is empty");
    initialize(T0);
}

template <class T>
Gsr_t<T>::Gsr_t(const Handle<YieldTermStructure_t<T> > &termStructure,
                const std::vector<Date> &volstepdates,
                const std::vector<Handle<Quote_t<T> > > &volatilities,
                const std::vector<Handle<Quote_t<T> > > &reversions,
                const Time T0)
    : Gaussian1dModel(termStructure), CalibratedModel_t<T>(2),
      reversion_(this->arguments_[0]), sigma_(this->arguments_[1]),
      volatilities_(volatilities), reversions_(reversions),
      volstepdates_(volstepdates) {

    QL_REQUIRE(!termStructure.empty(), "yield term structure handle is empty");
    initialize(T0);
}

template <class T> void Gsr_t<T>::updateTimes() const {
    volsteptimes_.clear();
    int j = 0;
    for (std::vector<Date>::const_iterator i = volstepdates_.begin();
         i != volstepdates_.end(); ++i, ++j) {
        volsteptimes_.push_back(this->termStructure()->timeFromReference(*i));
        volsteptimesArray_[j] = volsteptimes_[j];
        if (j == 0)
            QL_REQUIRE(volsteptimes_[0] > 0.0, "volsteptimes must be positive ("
                                                   << volsteptimes_[0] << ")");
        else
            QL_REQUIRE(volsteptimes_[j] > volsteptimes_[j - 1],
                       "volsteptimes must be strictly increasing ("
                           << volsteptimes_[j - 1] << "@" << (j - 1) << ", "
                           << volsteptimes_[j] << "@" << j << ")");
    }
    if (this->stateProcess_ != NULL)
        boost::static_pointer_cast<GsrProcess_t<T> >(this->stateProcess_)
            ->flushCache();
}

template <class T> void Gsr_t<T>::updateState() const {
    for (Size i = 0; i < sigma_.size(); i++) {
        sigma_.setParam(i, volatilities_[i]->value());
    }
    for (Size i = 0; i < reversion_.size(); i++) {
        reversion_.setParam(i, reversions_[i]->value());
    }
    boost::static_pointer_cast<GsrProcess_t<T> >(this->stateProcess_)
        ->flushCache();
}

template <class T> void Gsr_t<T>::initialize(Time T0) {

    volsteptimesArray_ = Array_t<Time>(volstepdates_.size());

    updateTimes();

    QL_REQUIRE(volatilities_.size() == volsteptimes_.size() + 1,
               "there must be n+1 volatilities ("
                   << volatilities_.size() << ") for n volatility step times ("
                   << volsteptimes_.size() << ")");
    // sigma_ =
    // PiecewiseConstantParameter(volsteptimes_,PositiveConstraint());
    sigma_ =
        PiecewiseConstantParameter_t<T>(volsteptimes_, NoConstraint_t<T>());

    QL_REQUIRE(reversions_.size() == 1 ||
                   reversions_.size() == volsteptimes_.size() + 1,
               "there must be 1 or n+1 reversions ("
                   << reversions_.size() << ") for n volatility step times ("
                   << volsteptimes_.size() << ")");
    if (reversions_.size() == 1) {
        reversion_ = ConstantParameter_t<T>(reversions_[0]->value(),
                                            NoConstraint_t<T>());
    } else {
        reversion_ =
            PiecewiseConstantParameter_t<T>(volsteptimes_, NoConstraint_t<T>());
    }

    this->stateProcess_ =
        boost::shared_ptr<GsrProcess_t<T> >(new GsrProcess_t<T>(
            volsteptimesArray_, sigma_.params(), reversion_.params(), T0));

    this->registerWith(this->termStructure());

    this->registerWith(this->stateProcess_);
    for (Size i = 0; i < reversions_.size(); ++i)
        this->registerWith(reversions_[i]);

    for (Size i = 0; i < volatilities_.size(); ++i)
        this->registerWith(volatilities_[i]);

    updateState();
}

template <class T>
const T
Gsr_t<T>::zerobondImpl(const Time T0, const Time t, const T y,
                       const Handle<YieldTermStructure_t<T> > &yts) const {

    this->calculate();

    if (t == 0.0)
        return yts.empty() ? this->termStructure()->discount(T0, true)
                           : yts->discount(T0, true);

    boost::shared_ptr<GsrProcess_t<T> > p =
        boost::dynamic_pointer_cast<GsrProcess_t<T> >(this->stateProcess_);

    T x = y * p->stdDeviation(0.0, 0.0, t) +
          this->stateProcess_->expectation(0.0, 0.0, t);
    T gtT = p->G(t, T0, x);

    T d = yts.empty()
              ? this->termStructure()->discount(T0, true) /
                    this->termStructure()->discount(t, true)
              : yts->discount(T0, true) / yts->discount(t, true);

    return d * QLFCT::exp(-x * gtT - 0.5 * p->y(t) * gtT * gtT);
}

template <class T>
const T
Gsr_t<T>::numeraireImpl(const Time t, const T y,
                        const Handle<YieldTermStructure_t<T> > &yts) const {

    this->calculate();

    boost::shared_ptr<GsrProcess_t<T> > p =
        boost::dynamic_pointer_cast<GsrProcess_t<T> >(this->stateProcess_);

    if (t == 0)
        return yts.empty()
                   ? this->termStructure()->discount(p->getForwardMeasureTime(),
                                                     true)
                   : yts->discount(p->getForwardMeasureTime());
    return this->zerobond(p->getForwardMeasureTime(), t, y, yts);
}

} // namespace QuantLib

#endif
