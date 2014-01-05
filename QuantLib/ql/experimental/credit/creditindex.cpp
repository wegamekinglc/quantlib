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

#include <ql/experimental/credit/creditindex.hpp>
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/settings.hpp>

namespace QuantLib {

    CreditIndex::CreditIndex(const std::string& familyName,
                             const Period& tenor,
                             Natural fixingDays,
                             DateGeneration::Rule dateRule, 
                             const Currency& currency,
                             const Calendar& fixingCalendar,
                             const DayCounter& dayCounter,
                             Seniority indexSeniority)
    : familyName_(familyName), 
      indexSeniority_(indexSeniority),
      tenor_(tenor), 
      fixingDays_(fixingDays), 
      fixingCalendar_(fixingCalendar),
      dayCounter_(dayCounter),
      dateRule_(dateRule) {
          registerWith(Settings::instance().evaluationDate());
          registerWith(IndexManager::instance().notifier(name()));
    }

    void CreditIndex::addFixing(const Date& fixingDate,
                                Real fixing,
                                bool forceOverwrite) {
        Index::addFixing(fixingDate, fixing, forceOverwrite);
    }

    Date CreditIndex::valueDate(const Date& fixingDate) const {
        QL_REQUIRE(isValidFixingDate(fixingDate),
                   "Fixing date " << fixingDate << " is not valid");
        return fixingCalendar().advance(fixingDate, fixingDays_, Days);
    }

    Date CreditIndex::maturityDate(const Date& valDate) const {
        return fixingCalendar().advance(valDate, tenor_, Unadjusted);
    }

    // Single Name index definitions -----------------------------------------

    Rate SingleNameCreditIndex::forecastFixing(const Date& fixingDate, 
        Real recovery) const {
        const boost::shared_ptr<CreditDefaultSwap> cdsSwap =
            underlyingSwap(fixingDate);
        cdsSwap->setPricingEngine(boost::shared_ptr<PricingEngine>(
            new MidPointCdsEngine(defaultProbTermStructure(), recovery, 
                discountTermStructure_)));
        return cdsSwap->fairSpread();
    }

    /* A single name index will return a fixing even if there are default 
      events triggering its contractual conditions. The knockability is  
      taken to be a concept associated to the coupons only.
    */
    Rate SingleNameCreditIndex::fixing(
        const Date& fixingDate,
        bool forecastTodaysFixing) const {
        QL_REQUIRE(isValidFixingDate(fixingDate),
                   "Fixing date " << fixingDate << " is not valid");
        Date today = Settings::instance().evaluationDate();
        // keeps the IR logic
        bool enforceTodaysHistoricFixings =
            Settings::instance().enforcesTodaysHistoricFixings();
        if(fixingDate < today ||
            (fixingDate == today && enforceTodaysHistoricFixings 
             && !forecastTodaysFixing)) { 
            Real fixingVal =
                IndexManager::instance().getHistory(name())[fixingDate];
            QL_REQUIRE(fixingVal != Null<Real>(),
                "Missing " << name() << " fixing for " << fixingDate);
            return fixingVal;
        }else if(fixingDate == today && !forecastTodaysFixing){
            // might have been fixed
            try {
                Rate pastFixing =
                    IndexManager::instance().getHistory(name())[fixingDate];
                if (pastFixing != Null<Real>())
                    return pastFixing;
                else
                    ;   // fall through and forecast
            } catch (Error&) {
                ;       // fall through and forecast
            }
        }
        return forecastFixing(fixingDate);
    }

    boost::shared_ptr<CreditDefaultSwap> 
        SingleNameCreditIndex::underlyingSwap(const Date& fixingDate) const {
            Schedule cdsSchedule =
                MakeSchedule()
                    .from(fixingDate)
                    .to(fixingDate+tenor_)
                    .withFrequency(paymentFreq_)
                    .withCalendar(fixingCalendar_)
                    .withTerminationDateConvention(Unadjusted)
                    .withRule(dateRule_);
            boost::shared_ptr<CreditDefaultSwap> cdsSwap(new
                CreditDefaultSwap(Protection::Buyer, 1., 1., cdsSchedule,
                    Following, dayCounter_));
            return cdsSwap;
    }

};
