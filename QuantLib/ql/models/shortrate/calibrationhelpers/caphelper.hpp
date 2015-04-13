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

/*! \file caphelper.hpp
    \brief CapHelper calibration helper
*/

#ifndef quantlib_cap_calibration_helper_hpp
#define quantlib_cap_calibration_helper_hpp

#include <ql/models/calibrationhelper.hpp>
#include <ql/instruments/capfloor.hpp>
#include <ql/pricingengines/capfloor/blackcapfloorengine.hpp>
#include <ql/pricingengines/capfloor/discretizedcapfloor.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/time/schedule.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/cashflows/cashflowvectors.hpp>

namespace QuantLib {

//! calibration helper for ATM cap

template <class T> class CapHelper_t : public CalibrationHelper_t<T> {
  public:
    CapHelper_t(const Period &length, const Handle<Quote_t<T> > &volatility,
                const boost::shared_ptr<IborIndex> &index,
                // data for ATM swap-rate calculation
                Frequency fixedLegFrequency,
                const DayCounter &fixedLegDayCounter, bool includeFirstSwaplet,
                const Handle<YieldTermStructure_t<T> > &termStructure,
                typename CalibrationHelper_t<T>::CalibrationErrorType
                    errorType = CalibrationHelper_t<T>::RelativePriceError);
    virtual void addTimesTo(std::list<Time> &times) const;
    virtual T modelValue() const;
    virtual T blackPrice(T volatility) const;

  private:
    void performCalculations() const;
    mutable boost::shared_ptr<Cap> cap_;
    const Period length_;
    const boost::shared_ptr<IborIndex> index_;
    const Frequency fixedLegFrequency_;
    const DayCounter fixedLegDayCounter_;
    const bool includeFirstSwaplet_;
};

typedef CapHelper_t<Real> CapHelper;

// implementation

template <class T>
CapHelper_t<T>::CapHelper_t(
    const Period &length, const Handle<Quote_t<T> > &volatility,
    const boost::shared_ptr<IborIndex> &index, Frequency fixedLegFrequency,
    const DayCounter &fixedLegDayCounter, bool includeFirstSwaplet,
    const Handle<YieldTermStructure_t<T> > &termStructure,
    typename CalibrationHelper_t<T>::CalibrationErrorType errorType)
    : CalibrationHelper_t<T>(volatility, termStructure, errorType),
      length_(length), index_(index), fixedLegFrequency_(fixedLegFrequency),
      fixedLegDayCounter_(fixedLegDayCounter),
      includeFirstSwaplet_(includeFirstSwaplet) {

    this->registerWith(index_);
}

template <class T>
void CapHelper_t<T>::addTimesTo(std::list<Time> &times) const {
    this->calculate();
    typename CapFloor_t<T>::arguments args;
    cap_->setupArguments(&args);
    std::vector<Time> capTimes =
        DiscretizedCapFloor(args, this->termStructure_->referenceDate(),
                            this->termStructure_->dayCounter())
            .mandatoryTimes();
    times.insert(times.end(), capTimes.begin(), capTimes.end());
}

template <class T> T CapHelper_t<T>::modelValue() const {
    this->calculate();
    cap_->setPricingEngine(this->engine_);
    return cap_->NPV();
}

template <class T> T CapHelper_t<T>::blackPrice(T sigma) const {
    this->calculate();
    boost::shared_ptr<Quote_t<T> > vol(new SimpleQuote_t<T>(sigma));
    boost::shared_ptr<PricingEngine> black(new BlackCapFloorEngine_t<T>(
        this->termStructure_, Handle<Quote_t<T> >(vol)));
    cap_->setPricingEngine(black);
    T value = cap_->NPV();
    cap_->setPricingEngine(this->engine_);
    return value;
}

template <class T> void CapHelper_t<T>::performCalculations() const {

    Period indexTenor = index_->tenor();
    T fixedRate = 0.04; // dummy value
    Date startDate, maturity;
    if (includeFirstSwaplet_) {
        startDate = this->termStructure_->referenceDate();
        maturity = this->termStructure_->referenceDate() + length_;
    } else {
        startDate = this->termStructure_->referenceDate() + indexTenor;
        maturity = this->termStructure_->referenceDate() + length_;
    }
    boost::shared_ptr<IborIndex> dummyIndex(new IborIndex(
        "dummy", indexTenor, index_->fixingDays(), index_->currency(),
        index_->fixingCalendar(), index_->businessDayConvention(),
        index_->endOfMonth(), this->termStructure_->dayCounter(),
        this->termStructure_));

    std::vector<T> nominals(1, 1.0);

    Schedule floatSchedule(
        startDate, maturity, index_->tenor(), index_->fixingCalendar(),
        index_->businessDayConvention(), index_->businessDayConvention(),
        DateGeneration::Forward, false);
    Leg floatingLeg =
        IborLeg(floatSchedule, index_)
            .withNotionals(nominals)
            .withPaymentAdjustment(index_->businessDayConvention())
            .withFixingDays(0);

    Schedule fixedSchedule(startDate, maturity, Period(fixedLegFrequency_),
                           index_->fixingCalendar(), Unadjusted, Unadjusted,
                           DateGeneration::Forward, false);
    Leg fixedLeg = FixedRateLeg(fixedSchedule)
                       .withNotionals(nominals)
                       .withCouponRates(fixedRate, fixedLegDayCounter_)
                       .withPaymentAdjustment(index_->businessDayConvention());

    Swap swap(floatingLeg, fixedLeg);
    swap.setPricingEngine(boost::shared_ptr<PricingEngine>(
        new DiscountingSwapEngine_t<T>(this->termStructure_, false)));
    T fairRate = fixedRate - swap.NPV() / (swap.legBPS(1) / 1.0e-4);
    cap_ = boost::shared_ptr<Cap_t<T> >(
        new Cap_t<T>(floatingLeg, std::vector<T>(1, fairRate)));

    CalibrationHelper_t<T>::performCalculations();
}
}

#endif
