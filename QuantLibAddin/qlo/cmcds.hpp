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

#ifndef qla_cmcds_hpp
#define qla_cmcds_hpp

#include <qlo/index.hpp>
#include <qlo/schedule.hpp>
#include <qlo/leg.hpp>
#include <qlo/baseinstruments.hpp>
#include <qlo/pricingengines.hpp>

#include <ql/time/schedule.hpp>
#include <ql/time/dategenerationrule.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/default.hpp>

#include <ql/experimental/credit/issuer.hpp>
#include <ql/experimental/credit/creditindex.hpp>
#include <ql/experimental/credit/creditcmcoupon.hpp>
#include <ql/experimental/credit/creditcouponpricer.hpp>

namespace QuantLibAddin {

    class Quote;
    class RecoveryRateQuote;

    OH_OBJ_CLASS(CreditIndex, Index);

    class SingleNameCreditIndex : public CreditIndex {
    public:
        SingleNameCreditIndex(
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

              bool permanent);
      protected:
        OH_OBJ_CTOR(SingleNameCreditIndex, CreditIndex);
    };


    class CmCdsLeg : public Leg {
      public:
        CmCdsLeg(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            QuantLib::BusinessDayConvention paymentConvention,
            const std::vector<QuantLib::Real>& nominals,
            const boost::shared_ptr<QuantLib::Schedule>& schedule,
            const std::vector<QuantLib::Natural>& fixingDays,
            bool isInArrears,
            const QuantLib::DayCounter& paymentDayCounter,
            //const std::vector<QuantLib::Rate>& floors,
            const std::vector<QuantLib::Real>& gearings,
            const boost::shared_ptr<QuantLib::SingleNameCreditIndex>& index,
            const std::vector<QuantLib::Spread>& spreads,
            //const std::vector<QuantLib::Rate>& caps,
            bool permanent);
    };


    OH_LIB_CLASS(CreditCouponPricer, QuantLib::CreditCouponPricer);


    class CdsCmCouponPricer : public CreditCouponPricer {
    public:
        CdsCmCouponPricer(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const QuantLib::Handle<QuantLib::Quote>& vol,
            const QuantLib::Handle<QuantLib::RecoveryRateQuote>& recovery,
            bool permanent);
    protected:
        OH_OBJ_CTOR(CdsCmCouponPricer, CreditCouponPricer);
    };


    class ConstantMaturityCDS : public Instrument {
    public:
        ConstantMaturityCDS(
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
            bool permanent);

    };


    class BalckConstantMaturityCDSEngine: public PricingEngine {
    public:
        BalckConstantMaturityCDSEngine(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>&,
            QuantLib::Real recoveryRate,
            const QuantLib::Handle<QuantLib::YieldTermStructure>&,
            const QuantLib::Handle<QuantLib::Quote>& vol,
            bool permanent);
    };

}

#endif
