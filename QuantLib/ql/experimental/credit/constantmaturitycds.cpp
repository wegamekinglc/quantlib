/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013 Jose Aparicio

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

#include <ql/experimental/credit/constantmaturitycds.hpp>
#include <ql/instruments/claim.hpp>


namespace QuantLib {

    ConstantMaturityCDS::ConstantMaturityCDS(
        Protection::Side side,
        Real notional,
        const Schedule& schedule,
        const boost::shared_ptr<SingleNameCreditIndex>& creditIndex,
        BusinessDayConvention convention,
        const DayCounter& dayCounter,
        Natural fixingDays,
        Real gearing,// rate , in [0,1]
        Rate cap, // only one 
        bool settlesAccrual,
        bool settlesAtDefT,
        const Date& protectionStart,
        const boost::shared_ptr<Claim>& claim
        )
    : side_(side), notional_(notional), creditIndex_(creditIndex), 
      paymentConvention_(convention), dayCount_(dayCounter), 
      gearing_(gearing), settlesAccrual_(settlesAccrual), 
      paysAtDefaultTime_(settlesAtDefT), 
      claim_(claim),
      protectionStart_(protectionStart == Null<Date>() ? schedule[0] :
                                                         protectionStart)
    {
        // protection start should fall within the first coupon period, 
        //   it is really an inception/settlement date....
        QL_REQUIRE(gearing>=0. , "Incorrect gearing value.");

        // create cmcds coupon leg
        couponLeg_ = CmCdsLeg(schedule, creditIndex)
            .withNotionals(notional)
            .withPaymentDayCounter(dayCounter)
            .withPaymentAdjustment(convention)
            .withFixingDays(fixingDays)
            .withGearings(gearing)
            .withCaps(cap)
            .withAccrualSettlement(settlesAccrual)
			.withGlobalProtectionDate(protectionStart_);

        if (!claim_)
            claim_ = boost::shared_ptr<Claim>(new FaceValueClaim);
        registerWith(claim_);

    }


    bool ConstantMaturityCDS::isExpired() const {
        for (Leg::const_reverse_iterator i = couponLeg_.rbegin();
                                         i != couponLeg_.rend(); ++i) {
            if (!(*i)->hasOccurred())
                return false;
        }
        return true;
    }

    void ConstantMaturityCDS::setupExpired() const {
        Instrument::setupExpired();
    }

    ConstantMaturityCDS::arguments::arguments()
    : side(Protection::Side(-1)), notional(Null<Real>()),
      gearing(Null<Rate>()), protectionStart(Date()) {}


    void ConstantMaturityCDS::setupArguments(
                                       PricingEngine::arguments* args) const {
        ConstantMaturityCDS::arguments* arguments =
            dynamic_cast<ConstantMaturityCDS::arguments*>(args);
        QL_REQUIRE(arguments != 0, "wrong argument type");

        arguments->side = side_;
        arguments->notional = notional_;
        arguments->leg = couponLeg_;
        arguments->settlesAccrual = settlesAccrual_;
        arguments->paysAtDefaultTime = paysAtDefaultTime_;
        arguments->claim = claim_;
        arguments->gearing = gearing_;
        arguments->protectionStart = protectionStart_;
        arguments->creditIndex = creditIndex_;
    }

    void ConstantMaturityCDS::arguments::validate() const {
        QL_REQUIRE(side != Protection::Side(-1), "side not set");
        QL_REQUIRE(notional != Null<Real>(), "notional not set");
        QL_REQUIRE(notional != 0.0, "null notional set");
        QL_REQUIRE(!leg.empty(), "coupons not set");
        QL_REQUIRE(claim, "claim not set");
        QL_REQUIRE(creditIndex, "index not set");
        QL_REQUIRE(protectionStart != Null<Date>(),
                   "protection start date not set");
        QL_REQUIRE(gearing != Null<Real>(), "gearing not set");
    }

    void ConstantMaturityCDS::fetchResults(
                                      const PricingEngine::results* r) const {
        Instrument::fetchResults(r);

        const ConstantMaturityCDS::results* results =
            dynamic_cast<const ConstantMaturityCDS::results*>(r);
        QL_REQUIRE(results != 0, "wrong result type");

        fairGearing_ = results->fairGearingFactor;
        couponLegNPV_ = results->couponLegNPV;
        defaultLegNPV_ = results->defaultLegNPV;
    }

    void ConstantMaturityCDS::results::reset() {
        Instrument::results::reset();
        couponLegNPV = Null<Real>();
        defaultLegNPV = Null<Real>();
        fairGearingFactor = Null<Real>();
        value = Null<Real>();
    }

 
    Real ConstantMaturityCDS::couponLegNPV() const {
        calculate();
        QL_REQUIRE(couponLegNPV_ != Null<Rate>(),
                   "coupon-leg NPV not available");
        return couponLegNPV_;
    }

    Real ConstantMaturityCDS::defaultLegNPV() const {
        calculate();
        QL_REQUIRE(defaultLegNPV_ != Null<Rate>(),
                   "default-leg NPV not available");
        return defaultLegNPV_;
    }

    Real ConstantMaturityCDS::fairGearing() const {
        calculate();
        QL_REQUIRE(fairGearing_ != Null<Real>(),
                   "fair gearing not available");
        return fairGearing_;
    }

}
