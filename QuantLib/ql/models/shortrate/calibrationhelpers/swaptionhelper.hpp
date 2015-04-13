/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2007 StatPro Italia srl
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

/*! \file swaptionhelper.hpp
    \brief Swaption calibration helper
*/

#ifndef quantlib_swaption_calibration_helper_hpp
#define quantlib_swaption_calibration_helper_hpp

#include <ql/models/calibrationhelper.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>
#include <ql/pricingengines/swaption/discretizedswaption.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/time/schedule.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/indexes/iborindex.hpp>

namespace QuantLib {

//! calibration helper for ATM swaption

template <class T> class SwaptionHelper_t : public CalibrationHelper_t<T> {
  public:
    SwaptionHelper_t(const Period &maturity, const Period &length,
                     const Handle<Quote_t<T> > &volatility,
                     const boost::shared_ptr<IborIndex_t<T> > &index,
                     const Period &fixedLegTenor,
                     const DayCounter &fixedLegDayCounter,
                     const DayCounter &floatingLegDayCounter,
                     const Handle<YieldTermStructure_t<T> > &termStructure,
                     typename CalibrationHelper_t<T>::CalibrationErrorType
                         errorType = CalibrationHelper_t<T>::RelativePriceError,
                     const T strike = Null<T>(), const T nominal = 1.0);

    SwaptionHelper_t(const Date &exerciseDate, const Period &length,
                     const Handle<Quote_t<T> > &volatility,
                     const boost::shared_ptr<IborIndex_t<T> > &index,
                     const Period &fixedLegTenor,
                     const DayCounter &fixedLegDayCounter,
                     const DayCounter &floatingLegDayCounter,
                     const Handle<YieldTermStructure_t<T> > &termStructure,
                     typename CalibrationHelper_t<T>::CalibrationErrorType
                         errorType = CalibrationHelper_t<T>::RelativePriceError,
                     const T strike = Null<T>(), const T nominal = 1.0);

    SwaptionHelper_t(const Date &exerciseDate, const Date &endDate,
                     const Handle<Quote_t<T> > &volatility,
                     const boost::shared_ptr<IborIndex_t<T> > &index,
                     const Period &fixedLegTenor,
                     const DayCounter &fixedLegDayCounter,
                     const DayCounter &floatingLegDayCounter,
                     const Handle<YieldTermStructure_t<T> > &termStructure,
                     typename CalibrationHelper_t<T>::CalibrationErrorType
                         errorType = CalibrationHelper_t<T>::RelativePriceError,
                     const T strike = Null<T>(), const T nominal = 1.0);

    virtual void addTimesTo(std::list<Time> &times) const;
    virtual T modelValue() const;
    virtual T blackPrice(T volatility) const;

    boost::shared_ptr<VanillaSwap_t<T> > underlyingSwap() const {
        this->calculate();
        return swap_;
    }
    boost::shared_ptr<Swaption_t<T> > swaption() const {
        this->calculate();
        return swaption_;
    }

  private:
    void performCalculations() const;
    mutable Date exerciseDate_, endDate_;
    const Period maturity_, length_, fixedLegTenor_;
    const boost::shared_ptr<IborIndex_t<T> > index_;
    const DayCounter fixedLegDayCounter_, floatingLegDayCounter_;
    const T strike_, nominal_;
    mutable T exerciseRate_;
    mutable boost::shared_ptr<VanillaSwap_t<T> > swap_;
    mutable boost::shared_ptr<Swaption_t<T> > swaption_;
};

typedef SwaptionHelper_t<Real> SwaptionHelper;

// implementation

template <class T>
SwaptionHelper_t<T>::SwaptionHelper_t(
    const Period &maturity, const Period &length,
    const Handle<Quote_t<T> > &volatility,
    const boost::shared_ptr<IborIndex_t<T> > &index,
    const Period &fixedLegTenor, const DayCounter &fixedLegDayCounter,
    const DayCounter &floatingLegDayCounter,
    const Handle<YieldTermStructure_t<T> > &termStructure,
    typename CalibrationHelper_t<T>::CalibrationErrorType errorType,
    const T strike, const T nominal)

    : CalibrationHelper_t<T>(volatility, termStructure, errorType),
      exerciseDate_(Null<Date>()), endDate_(Null<Date>()), maturity_(maturity),
      length_(length), fixedLegTenor_(fixedLegTenor), index_(index),
      fixedLegDayCounter_(fixedLegDayCounter),
      floatingLegDayCounter_(floatingLegDayCounter), strike_(strike),
      nominal_(nominal) {

    this->registerWith(index_);
}

template <class T>
SwaptionHelper_t<T>::SwaptionHelper_t(
    const Date &exerciseDate, const Period &length,
    const Handle<Quote_t<T> > &volatility,
    const boost::shared_ptr<IborIndex_t<T> > &index,
    const Period &fixedLegTenor, const DayCounter &fixedLegDayCounter,
    const DayCounter &floatingLegDayCounter,
    const Handle<YieldTermStructure_t<T> > &termStructure,
    typename CalibrationHelper_t<T>::CalibrationErrorType errorType,
    const T strike, const T nominal)
    : CalibrationHelper_t<T>(volatility, termStructure, errorType),
      exerciseDate_(exerciseDate), endDate_(Null<Date>()), maturity_(0 * Days),
      length_(length), fixedLegTenor_(fixedLegTenor), index_(index),
      fixedLegDayCounter_(fixedLegDayCounter),
      floatingLegDayCounter_(floatingLegDayCounter), strike_(strike),
      nominal_(nominal) {

    this->registerWith(index_);
}

template <class T>
SwaptionHelper_t<T>::SwaptionHelper_t(
    const Date &exerciseDate, const Date &endDate,
    const Handle<Quote_t<T> > &volatility,
    const boost::shared_ptr<IborIndex_t<T> > &index,
    const Period &fixedLegTenor, const DayCounter &fixedLegDayCounter,
    const DayCounter &floatingLegDayCounter,
    const Handle<YieldTermStructure_t<T> > &termStructure,
    typename CalibrationHelper_t<T>::CalibrationErrorType errorType,
    const T strike, const T nominal)
    : CalibrationHelper_t<T>(volatility, termStructure, errorType),
      exerciseDate_(exerciseDate), endDate_(endDate), maturity_(0 * Days),
      length_(0 * Days), fixedLegTenor_(fixedLegTenor), index_(index),
      fixedLegDayCounter_(fixedLegDayCounter),
      floatingLegDayCounter_(floatingLegDayCounter), strike_(strike),
      nominal_(nominal) {

    this->registerWith(index_);
}

template <class T>
void SwaptionHelper_t<T>::addTimesTo(std::list<Time> &times) const {
    this->calculate();
    typename Swaption_t<T>::arguments args;
    swaption_->setupArguments(&args);
    std::vector<Time> swaptionTimes =
        DiscretizedSwaption_t<T>(args, this->termStructure_->referenceDate(),
                            this->termStructure_->dayCounter())
            .mandatoryTimes();
    times.insert(times.end(), swaptionTimes.begin(), swaptionTimes.end());
}

template <class T> T SwaptionHelper_t<T>::modelValue() const {
    this->calculate();
    swaption_->setPricingEngine(this->engine_);
    return swaption_->NPV();
}

template <class T> T SwaptionHelper_t<T>::blackPrice(T sigma) const {
    this->calculate();
    Handle<Quote_t<T> > vol(
        boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(sigma)));
    boost::shared_ptr<PricingEngine> black(
        new BlackSwaptionEngine_t<T>(this->termStructure_, vol));
    swaption_->setPricingEngine(black);
    T value = swaption_->NPV();
    swaption_->setPricingEngine(this->engine_);
    return value;
}

template <class T> void SwaptionHelper_t<T>::performCalculations() const {

    Calendar calendar = index_->fixingCalendar();
    Natural fixingDays = index_->fixingDays();

    Date exerciseDate = exerciseDate_;
    if (exerciseDate == Null<Date>())
        exerciseDate =
            calendar.advance(this->termStructure_->referenceDate(), maturity_,
                             index_->businessDayConvention());

    Date startDate = calendar.advance(exerciseDate, fixingDays, Days,
                                      index_->businessDayConvention());

    Date endDate = endDate_;
    if (endDate == Null<Date>())
        endDate = calendar.advance(startDate, length_,
                                   index_->businessDayConvention());

    Schedule fixedSchedule(startDate, endDate, fixedLegTenor_, calendar,
                           index_->businessDayConvention(),
                           index_->businessDayConvention(),
                           DateGeneration::Forward, false);
    Schedule floatSchedule(startDate, endDate, index_->tenor(), calendar,
                           index_->businessDayConvention(),
                           index_->businessDayConvention(),
                           DateGeneration::Forward, false);

    boost::shared_ptr<PricingEngine> swapEngine(
        new DiscountingSwapEngine_t<T>(this->termStructure_, false));

    typename VanillaSwap_t<T>::Type type = VanillaSwap_t<T>::Receiver;

    VanillaSwap_t<T> temp(VanillaSwap_t<T>::Receiver, nominal_, fixedSchedule,
                          0.0, fixedLegDayCounter_, floatSchedule, index_, 0.0,
                          floatingLegDayCounter_);
    temp.setPricingEngine(swapEngine);
    T forward = temp.fairRate();
    if (strike_ == Null<T>()) {
        exerciseRate_ = forward;
    } else {
        exerciseRate_ = strike_;
        type = strike_ <= forward ? VanillaSwap_t<T>::Receiver
                                  : VanillaSwap_t<T>::Payer;
        // ensure that calibration instrument is out of the money
    }
    swap_ = boost::shared_ptr<VanillaSwap_t<T> >(new VanillaSwap_t<T>(
        type, nominal_, fixedSchedule, exerciseRate_, fixedLegDayCounter_,
        floatSchedule, index_, 0.0, floatingLegDayCounter_));
    swap_->setPricingEngine(swapEngine);

    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    swaption_ =
        boost::shared_ptr<Swaption_t<T> >(new Swaption_t<T>(swap_, exercise));

    CalibrationHelper_t<T>::performCalculations();
}
}

#endif
