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

#include <qlo/cmcds.hpp>

#include <ql/currencies/europe.hpp>
#include <ql/experimental/credit/constantmaturitycds.hpp>
#include <ql/experimental/credit/blackcmcdsengine.hpp>

namespace QuantLibAddin {

    SingleNameCreditIndex::SingleNameCreditIndex(
              const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
              const std::string& familyName,
              // turn to a Handle if Issuer becomes observable
              const boost::shared_ptr<QuantLib::Issuer>& issuer, 
       //////       const QuantLib::DefaultProbKey& defKey, // standard
              const QuantLib::Period& tenor, 
              const QuantLib::Frequency payFreq,  // Quarterly
              QuantLib::Natural fixingDays, // 0 
              QuantLib::DateGeneration::Rule dateRule, // CDS
              const QuantLib::Calendar& fixingCalendar,
              const QuantLib::DayCounter& dayCounter,
              const QuantLib::Handle<QuantLib::YieldTermStructure>& yts,

              const std::vector<QuantLib::Date>& fixingDates,
              const std::vector<QuantLib::Real>& fixingRates,

              bool permanent
              )
    : CreditIndex(properties, permanent) 
    {
        const QuantLib::NorthAmericaCorpDefaultKey 
            defKey(QuantLib::EURCurrency(),
                   QuantLib::SeniorSec, 
                   QuantLib::Period(),
                   1. // amount threshold
                   );
        libraryObject_ = boost::shared_ptr<QuantLib::SingleNameCreditIndex>(
            new QuantLib::SingleNameCreditIndex(
                familyName, *issuer, defKey, tenor, payFreq, 
                fixingDays, dateRule, fixingCalendar, dayCounter, yts));
        // add fixings if any
        QL_REQUIRE(fixingDates.size() == fixingRates.size(), 
            "Fixing dates and rates arrays have different sizes.");
        for(QuantLib::Size i=0; i<fixingDates.size(); i++)
            libraryObject_->addFixing(fixingDates[i], fixingRates[i]);

    }


    CmCdsLeg::CmCdsLeg(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            QuantLib::BusinessDayConvention paymentConvention,
            const std::vector<QuantLib::Real>& nominals,
            const boost::shared_ptr<QuantLib::Schedule>& schedule,
            const std::vector<QuantLib::Natural>& fixingDays,
            bool isInArrears,
            const QuantLib::DayCounter& paymentDayCounter,
            const std::vector<QuantLib::Real>& gearings,
            const boost::shared_ptr<QuantLib::SingleNameCreditIndex>& index,
            const std::vector<QuantLib::Spread>& spreads,
            bool permanent)
    : Leg(properties, permanent)
    {
        libraryObject_ = boost::shared_ptr<QuantLib::Leg>(new
            QuantLib::Leg(QuantLib::CmCdsLeg(*schedule, index)
                .withNotionals(nominals)
                .withPaymentDayCounter(paymentDayCounter)
                .withPaymentAdjustment(paymentConvention)
                .withFixingDays(fixingDays)
                .withGearings(gearings)
                .withSpreads(spreads)
                ));
                
    }


    CdsCmCouponPricer::CdsCmCouponPricer(
        const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
        const QuantLib::Handle<QuantLib::Quote>& vol,
        const QuantLib::Handle<QuantLib::RecoveryRateQuote>& recovery,
        bool permanent)
    : CreditCouponPricer(properties, permanent) {
    libraryObject_ =
        boost::shared_ptr<QuantLib::CdsCmCouponPricer>(new
            QuantLib::CdsCmCouponPricer(vol, recovery));
    }


    ConstantMaturityCDS::ConstantMaturityCDS(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            QuantLib::Protection::Side side,
            QuantLib::Real notional,
            const boost::shared_ptr<QuantLib::Schedule>& schedule,
            QuantLib::Natural fixingDays,
            const boost::shared_ptr<QuantLib::SingleNameCreditIndex>& index,
            QuantLib::BusinessDayConvention paymentConvention,
            const QuantLib::DayCounter& dayCounter,
            QuantLib::Real gearing,
            QuantLib::Rate cap,
            bool settlesAccrual,
            bool paysAtDefT,
            const QuantLib::Date& protectionStart,
            bool permanent)
   : Instrument(properties, permanent) {
        libraryObject_ = boost::shared_ptr<QuantLib::ConstantMaturityCDS>(
            new QuantLib::ConstantMaturityCDS(
                side,
                notional,
                *schedule,
                index,
                paymentConvention,
                dayCounter,
                fixingDays, // fixing days
                gearing,
                cap,
                settlesAccrual,
                paysAtDefT,
                protectionStart,
                boost::shared_ptr<QuantLib::Claim>()));
    }


    BalckConstantMaturityCDSEngine::BalckConstantMaturityCDSEngine(
        const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
        const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& dfts,
        QuantLib::Real recoveryRate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& yts,
        const QuantLib::Handle<QuantLib::Quote>& vol,
        bool permanent)
    : PricingEngine(properties, permanent) {

        libraryObject_ = 
            boost::shared_ptr<QuantLib::BalckConstantMaturityCDSEngine>(
            new QuantLib::BalckConstantMaturityCDSEngine(dfts, recoveryRate, 
                                                         yts, vol));
    }


}
