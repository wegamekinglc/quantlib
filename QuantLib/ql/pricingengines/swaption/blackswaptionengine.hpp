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

/*! \file blackswaptionengine.hpp
    \brief Black-formula swaption engine
*/

#ifndef quantlib_pricers_black_swaption_hpp
#define quantlib_pricers_black_swaption_hpp

#include <ql/instruments/swaptionbase.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/quote.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/termstructures/volatility/swaption/swaptionconstantvol.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/exercise.hpp>

namespace QuantLib {

//! Black-formula swaption engine
/*! \ingroup swaptionengines

    \warning The engine assumes that the exercise date equals the
             start date of the passed swap.
*/
template <class T> class BlackSwaptionEngine_t : public Swaption_t<T>::engine {
  public:
    BlackSwaptionEngine_t(const Handle<YieldTermStructure_t<T> > &discountCurve,
                          T vol, const DayCounter &dc = Actual365Fixed(),
                          T displacement = 0.0);
    BlackSwaptionEngine_t(const Handle<YieldTermStructure_t<T> > &discountCurve,
                          const Handle<Quote_t<T> > &vol,
                          const DayCounter &dc = Actual365Fixed(),
                          T displacement = 0.0);
    BlackSwaptionEngine_t(const Handle<YieldTermStructure_t<T> > &discountCurve,
                          const Handle<SwaptionVolatilityStructure_t<T> > &vol,
                          T displacement = 0.0);
    void calculate() const;
    Handle<YieldTermStructure_t<T> > termStructure() { return discountCurve_; }
    Handle<SwaptionVolatilityStructure_t<T> > volatility() { return vol_; }
    T displacement() { return displacement_; }

  private:
    Handle<YieldTermStructure_t<T> > discountCurve_;
    Handle<SwaptionVolatilityStructure_t<T> > vol_;
    T displacement_;
};

typedef BlackSwaptionEngine_t<Real> BlackSwaptionEngine;

// implementation

template <class T>
BlackSwaptionEngine_t<T>::BlackSwaptionEngine_t(
    const Handle<YieldTermStructure_t<T> > &discountCurve, T vol,
    const DayCounter &dc, T displacement)
    : discountCurve_(discountCurve),
      vol_(boost::shared_ptr<SwaptionVolatilityStructure_t<T> >(
          new ConstantSwaptionVolatility_t<T>(0, NullCalendar(), Following, vol,
                                              dc))),
      displacement_(displacement) {
    this->registerWith(discountCurve_);
}

template <class T>
BlackSwaptionEngine_t<T>::BlackSwaptionEngine_t(
    const Handle<YieldTermStructure_t<T> > &discountCurve,
    const Handle<Quote_t<T> > &vol, const DayCounter &dc, T displacement)
    : discountCurve_(discountCurve),
      vol_(boost::shared_ptr<SwaptionVolatilityStructure_t<T> >(
          new ConstantSwaptionVolatility_t<T>(0, NullCalendar(), Following, vol,
                                              dc))),
      displacement_(displacement) {
    this->registerWith(discountCurve_);
    this->registerWith(vol_);
}

template <class T>
BlackSwaptionEngine_t<T>::BlackSwaptionEngine_t(
    const Handle<YieldTermStructure_t<T> > &discountCurve,
    const Handle<SwaptionVolatilityStructure_t<T> > &volatility, T displacement)
    : discountCurve_(discountCurve), vol_(volatility),
      displacement_(displacement) {
    this->registerWith(discountCurve_);
    this->registerWith(vol_);
}

template <class T> void BlackSwaptionEngine_t<T>::calculate() const {
    static const T basisPoint = 1.0e-4;

    Date exerciseDate = this->arguments_.exercise->date(0);

    // the part of the swap preceding exerciseDate should be truncated
    // to avoid taking into account unwanted cashflows
    VanillaSwap_t<T> swap = *this->arguments_.swap;

    T strike = swap.fixedRate();

    // using the discounting curve
    // swap.iborIndex() might be using a different forwarding curve
    swap.setPricingEngine(boost::shared_ptr<PricingEngine>(
        new DiscountingSwapEngine_t<T>(discountCurve_, false)));
    T atmForward = swap.fairRate();

    // Volatilities are quoted for zero-spreaded swaps.
    // Therefore, any spread on the floating leg must be removed
    // with a corresponding correction on the fixed leg.
    if (swap.spread() != 0.0) {
        T correction = swap.spread() *
                       QLFCT::abs(swap.floatingLegBPS() / swap.fixedLegBPS());
        strike -= correction;
        atmForward -= correction;
        this->results_.additionalResults["spreadCorrection"] = correction;
    } else {
        this->results_.additionalResults["spreadCorrection"] = 0.0;
    }
    this->results_.additionalResults["strike"] = strike;
    this->results_.additionalResults["atmForward"] = atmForward;

    // using the discounting curve
    swap.setPricingEngine(boost::shared_ptr<PricingEngine>(
        new DiscountingSwapEngine_t<T>(discountCurve_, false)));
    T annuity;
    switch (this->arguments_.settlementType) {
    case Settlement::Physical: {
        annuity = QLFCT::abs(swap.fixedLegBPS()) / basisPoint;
        break;
    }
    case Settlement::Cash: {
        const typename Leg_t<T>::Type &fixedLeg = swap.fixedLeg();
        boost::shared_ptr<FixedRateCoupon_t<T>> firstCoupon =
            boost::dynamic_pointer_cast<FixedRateCoupon_t<T>>(fixedLeg[0]);
        DayCounter dayCount = firstCoupon->dayCounter();
        T fixedLegCashBPS = CashFlows::bps(
            fixedLeg, InterestRate_t<T>(atmForward, dayCount, Compounded, Annual),
            false, discountCurve_->referenceDate());
        annuity = QLFCT::abs(fixedLegCashBPS / basisPoint);
        break;
    }
    default:
        QL_FAIL("unknown settlement type");
    }
    this->results_.additionalResults["annuity"] = annuity;

    // the swap length calculation might be improved using the value date
    // of the exercise date
    Time swapLength = vol_->swapLength(
        exerciseDate, this->arguments_.floatingPayDates.back());
    this->results_.additionalResults["swapLength"] = swapLength;

    T variance = vol_->blackVariance(exerciseDate, swapLength, strike);
    T stdDev = QLFCT::sqrt(variance);
    this->results_.additionalResults["stdDev"] = stdDev;
    Option::Type w = (this->arguments_.type == VanillaSwap_t<T>::Payer)
                         ? Option::Call
                         : Option::Put;
    this->results_.value =
        blackFormula(w, strike, atmForward, stdDev, annuity, displacement_);

    Time exerciseTime = vol_->timeFromReference(exerciseDate);
    this->results_.additionalResults["vega"] =
        QLFCT::sqrt(exerciseTime) *
        blackFormulaStdDevDerivative(strike, atmForward, stdDev, annuity,
                                     displacement_);
}
}

#endif
