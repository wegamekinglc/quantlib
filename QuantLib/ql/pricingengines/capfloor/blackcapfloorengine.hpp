/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 Ferdinando Ametrano
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2006 StatPro Italia srl
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

/*! \file blackcapfloorengine.hpp
    \brief Black-formula cap/floor engine
*/

#ifndef quantlib_pricers_black_capfloor_hpp
#define quantlib_pricers_black_capfloor_hpp

#include <ql/instruments/capfloorbase.hpp>
#include <ql/termstructures/volatility/optionlet/optionletvolatilitystructure.hpp>
#include <ql/quote.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/volatility/optionlet/constantoptionletvol.hpp>
#include <ql/time/calendars/nullcalendar.hpp>

namespace QuantLib {

//! Black-formula cap/floor engine
/*! \ingroup capfloorengines */
template <class T> class BlackCapFloorEngine_t : public CapFloor_t<T>::engine {
  public:
    BlackCapFloorEngine_t(const Handle<YieldTermStructure_t<T> > &discountCurve,
                          T vol, const DayCounter &dc = Actual365Fixed(),
                          T displacement = 0.0);
    BlackCapFloorEngine_t(const Handle<YieldTermStructure_t<T> > &discountCurve,
                          const Handle<Quote_t<T> > &vol,
                          const DayCounter &dc = Actual365Fixed(),
                          T displacement = 0.0);
    BlackCapFloorEngine_t(const Handle<YieldTermStructure_t<T> > &discountCurve,
                          const Handle<OptionletVolatilityStructure> &vol,
                          T displacement = 0.0);
    void calculate() const;
    Handle<YieldTermStructure_t<T> > termStructure() { return discountCurve_; }
    Handle<OptionletVolatilityStructure> volatility() { return vol_; }
    T displacement() { return displacement_; }

  private:
    Handle<YieldTermStructure_t<T> > discountCurve_;
    Handle<OptionletVolatilityStructure> vol_;
    T displacement_;
};

typedef BlackCapFloorEngine_t<Real> BlackCapFloorEngine;

// implementation

template <class T>
BlackCapFloorEngine_t<T>::BlackCapFloorEngine_t(
    const Handle<YieldTermStructure_t<T> > &discountCurve, T v,
    const DayCounter &dc, T displacement)
    : discountCurve_(discountCurve),
      vol_(boost::shared_ptr<OptionletVolatilityStructure>(
          new ConstantOptionletVolatility(0, NullCalendar(), Following, v,
                                          dc))),
      displacement_(displacement) {
    this->registerWith(discountCurve_);
}

template <class T>
BlackCapFloorEngine_t<T>::BlackCapFloorEngine_t(
    const Handle<YieldTermStructure_t<T> > &discountCurve,
    const Handle<Quote_t<T> > &v, const DayCounter &dc, T displacement)
    : discountCurve_(discountCurve),
      vol_(boost::shared_ptr<OptionletVolatilityStructure>(
          new ConstantOptionletVolatility(0, NullCalendar(), Following, v,
                                          dc))),
      displacement_(displacement) {
    this->registerWith(discountCurve_);
    this->registerWith(vol_);
}

template <class T>
BlackCapFloorEngine_t<T>::BlackCapFloorEngine_t(
    const Handle<YieldTermStructure_t<T> > &discountCurve,
    const Handle<OptionletVolatilityStructure> &volatility, T displacement)
    : discountCurve_(discountCurve), vol_(volatility),
      displacement_(displacement) {
    this->registerWith(discountCurve_);
    this->registerWith(vol_);
}

template <class T> void BlackCapFloorEngine_t<T>::calculate() const {
    T value = 0.0;
    T vega = 0.0;
    Size optionlets = this->arguments_.startDates.size();
    std::vector<T> values(optionlets, 0.0);
    std::vector<T> vegas(optionlets, 0.0);
    std::vector<T> stdDevs(optionlets, 0.0);
    typename CapFloor_t<T>::Type type = this->arguments_.type;
    Date today = vol_->referenceDate();
    Date settlement = discountCurve_->referenceDate();

    for (Size i = 0; i < optionlets; ++i) {
        Date paymentDate = this->arguments_.endDates[i];
        // handling of settlementDate, npvDate and includeSettlementFlows
        // should be implemented.
        // For the time being just discard expired caplets
        if (paymentDate > settlement) {
            DiscountFactor d = this->arguments_.nominals[i] * this->arguments_.gearings[i] *
                               discountCurve_->discount(paymentDate) *
                               this->arguments_.accrualTimes[i];

            T forward = this->arguments_.forwards[i];

            Date fixingDate = this->arguments_.fixingDates[i];
            Time sqrtTime = 0.0;
            if (fixingDate > today)
                sqrtTime = QLFCT::sqrt(vol_->timeFromReference(fixingDate));

            if (type == CapFloor_t<T>::Cap || type == CapFloor_t<T>::Collar) {
                T strike = this->arguments_.capRates[i];
                if (sqrtTime > 0.0) {
                    stdDevs[i] =
                        QLFCT::sqrt(vol_->blackVariance(fixingDate, strike));
                    vegas[i] = blackFormulaStdDevDerivative(strike, forward,
                                                            stdDevs[i], d,
                                                            displacement_) *
                               sqrtTime;
                }
                // include caplets with past fixing date
                values[i] = blackFormula(Option::Call, strike, forward,
                                         stdDevs[i], d, displacement_);
            }
            if (type == CapFloor_t<T>::Floor || type == CapFloor_t<T>::Collar) {
                T strike = this->arguments_.floorRates[i];
                T floorletVega = 0.0;
                if (sqrtTime > 0.0) {
                    stdDevs[i] =
                        QLFCT::sqrt(vol_->blackVariance(fixingDate, strike));
                    floorletVega = blackFormulaStdDevDerivative(strike, forward,
                                                                stdDevs[i], d,
                                                                displacement_) *
                                   sqrtTime;
                }
                T floorlet = blackFormula(Option::Put, strike, forward,
                                          stdDevs[i], d, displacement_);
                if (type == CapFloor_t<T>::Floor) {
                    values[i] = floorlet;
                    vegas[i] = floorletVega;
                } else {
                    // a collar is long a cap and short a floor
                    values[i] -= floorlet;
                    vegas[i] -= floorletVega;
                }
            }
            value += values[i];
            vega += vegas[i];
        }
    }
    this->results_.value = value;
    this->results_.additionalResults["vega"] = vega;

    this->results_.additionalResults["optionletsPrice"] = values;
    this->results_.additionalResults["optionletsVega"] = vegas;
    this->results_.additionalResults["optionletsAtmForward"] = this->arguments_.forwards;
    if (type != CapFloor_t<T>::Collar)
        this->results_.additionalResults["optionletsStdDev"] = stdDevs;
}
}

#endif
