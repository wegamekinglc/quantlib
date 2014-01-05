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

/*! \file creditcouponpricer.cpp
    \brief Pricing of a credit CM CDS.
*/

#include <ql/experimental/credit/creditcouponpricer.hpp>
#include <ql/experimental/credit/creditindex.hpp>
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/quotes/simplequote.hpp>

namespace QuantLib {

    CdsCmCouponPricer::CdsCmCouponPricer(Real vol,
                    Real recovery) 
    : vol_(boost::shared_ptr<Quote>(new SimpleQuote(vol))),
      recovery_(boost::shared_ptr<RecoveryRateQuote>(
        new RecoveryRateQuote(recovery))) { }

    Real CdsCmCouponPricer::swapletPrice() const {
        Date today = Settings::instance().evaluationDate();

        const Handle<DefaultProbabilityTermStructure>& defTS = 
            coupon_->creditIndex()->defaultProbTermStructure();
        const Handle<YieldTermStructure>& yTS =
            coupon_->creditIndex()->nominalTermStructure();

        // check for defaults; realized knock outs and realized accruals
        if(today >= refProtStart_) {
            boost::shared_ptr<DefaultEvent> defEvent = 
                coupon_->creditIndex()->issuer().defaultedBetween(
                refProtStart_, today, coupon_->creditIndex()->defaultKey());
            if(defEvent) {
                // Notice that if an accrual payment is to be generated now, 
                // the default lies within the accrual period and thence the  
                // coupon has fixed (thus no convexity adjustment)
                if(coupon_->settlesAccrual()
                    // event accrues only if in accrual period
                    && (coupon_->accrualStartDate() < defEvent->date()) ) {
                    return (gearing_ * effectiveRate_ + spread_) *
                        coupon_->accruedPeriod(defEvent->date()) * 
						yTS->discount(paymentDate_);
                }else{// knocked out without value
                    return 0. ;
                }
            }
        }

        // no defaults, compute expected accrual
        Real defaultAccrualTerm =0.;
        if(coupon_->settlesAccrual()) {
            Date effectiveAccrualDate = 
                coupon_->accrualStartDate() <= today ? 
                today : coupon_->accrualStartDate();
            // protection might start sooner but there would be no 
            // payments then.
            Date midPointDefaultDate = effectiveAccrualDate + 
                (coupon_->accrualEndDate() - effectiveAccrualDate) / 2 ;
			if(!coupon_->hasOccurred())
                defaultAccrualTerm = 
                /* in default within the subperiod leading to payments. 
				The coupon might be sensitive to former defaults but would 
				not generate the accrual payment*/
					defaultProb_ *
                // and no previous knock out during the protection period
                    defTS->survivalProbability(effectiveProtectStart_) * 
                    coupon_->dayCounter().yearFraction(
                        coupon_->accrualStartDate(), midPointDefaultDate,
                        coupon_->referencePeriodStart(),
	                    coupon_->referencePeriodEnd()) *
                    yTS->discount(paymentDate_);
        }

        Real convexAdjFactor = 1.;
        // Determine convexity adjustment:
        if(today < fixingDate_ ) {
            Real probsRatio = defTS->survivalProbability(refMaturity_) / 
                    defTS->survivalProbability(paymentDate_);
            convexAdjFactor = probsRatio + (1. - probsRatio) * 
                std::exp(std::pow(vol_->value(), 2.) * 
                         defTS->dayCounter().yearFraction(today, fixingDate_));
        }

		// expected value
        return 
            (survivalProb_ * accrualPeriod_ * discount_ + defaultAccrualTerm)*
            (gearing_ * effectiveRate_ * convexAdjFactor + spread_);
    }

    Rate CdsCmCouponPricer::capletPrice(const Rate effectiveCap) const {
        QL_REQUIRE(effectiveCap > 0., "Negative cap rate in CMCDS.");
        Date today = Settings::instance().evaluationDate();

        const Handle<DefaultProbabilityTermStructure>& defTS = 
            coupon_->creditIndex()->defaultProbTermStructure();
        const Handle<YieldTermStructure>& yTS =
            coupon_->creditIndex()->nominalTermStructure();

        Real capRate;

        // check for defaults:
        if(today >= refProtStart_) {
            boost::shared_ptr<DefaultEvent> defEvent = 
                coupon_->creditIndex()->issuer().defaultedBetween(
                refProtStart_, today, coupon_->creditIndex()->defaultKey());
            if(defEvent) {
                // Notice that if an accrual payment is to be generated now,
                // the default lies within the accrual period and thence the
                // coupon has fixed (thus no convexity adjustment)
                if(coupon_->settlesAccrual()
                    // event accrues only if in accrual period
                    && (coupon_->accrualStartDate() < defEvent->date()) ) {
                    return gearing_ * 
                        std::min(effectiveRate_ + spread_, effectiveCap)
						* coupon_->accruedPeriod(defEvent->date())
						* yTS->discount(paymentDate_);
                }else{// knocked out without value
                    return 0. ;
                }
            }
        }

        // no defaults, compute expected accrual
        Real defaultAccrualTerm =0.;
        if(coupon_->settlesAccrual()) {
            Date effectiveAccrualDate = 
                coupon_->accrualStartDate() <= today ? 
                today : coupon_->accrualStartDate();

            // protection might start sooner but there would be no 
            // payments then.
            Date midPointDefaultDate = effectiveAccrualDate + 
                (coupon_->accrualEndDate() - effectiveAccrualDate) / 2 ;
			if(!coupon_->hasOccurred())
                defaultAccrualTerm = 
                /* in default within the subperiod leading to payments. 
				The coupon might be sensitive to former defaults but would 
				not generate the accrual payment*/
					defaultProb_ *
                // and no previous knock out during the protection period
                    defTS->survivalProbability(effectiveProtectStart_) * 
                    coupon_->dayCounter().yearFraction(
                        coupon_->accrualStartDate(), midPointDefaultDate,
                        coupon_->referencePeriodStart(),
	                    coupon_->referencePeriodEnd()) *
                    yTS->discount(paymentDate_);
        }



        if(today < fixingDate_ ) {
            // reinventing the wheel? check reducing this to Blacks formula
            const Rate fwdRate = effectiveRate_;
            const Real vola = vol_->value(), vola2 = vola * vola;
            const Time t2fix = 
                defTS->dayCounter().yearFraction(today, fixingDate_);
            const Real phi1 = CumulativeNormalDistribution()(
                (std::log(effectiveCap/fwdRate)
                 - .5 * vola2*t2fix)
                / (vola * std::sqrt(t2fix)));
            const Real phi2 = CumulativeNormalDistribution()(
                (std::log(effectiveCap/fwdRate)
                 - 1.5 * vola2*t2fix)
                / (vola * std::sqrt(t2fix)));
            const Real phi3 = CumulativeNormalDistribution()(
                (std::log(fwdRate/effectiveCap)
                 - .5 * vola2*t2fix)
                / (vola * std::sqrt(t2fix)));
            const Real phi4 = CumulativeNormalDistribution()(
                (std::log(fwdRate/effectiveCap) 
                 + vola2*t2fix)
                / (vola * std::sqrt(t2fix)));
            Real probsRatio = defTS->survivalProbability(refMaturity_) / 
                    defTS->survivalProbability(paymentDate_);
            capRate = 
                fwdRate * probsRatio * phi1 + 
                fwdRate * (1.-probsRatio) * std::exp(vola2 * t2fix) * phi2 +
                effectiveCap * probsRatio * phi3 + 
                effectiveCap * (1.-probsRatio) * phi4;
            return (survivalProb_*discount_*coupon_->accrualPeriod() +
				defaultAccrualTerm) * gearing_ * capRate;
        }else{ // rate is fixed:
            capRate = std::min(effectiveCap, effectiveRate_);
            return (survivalProb_*discount_*coupon_->accrualPeriod() +
				defaultAccrualTerm) * gearing_ * capRate;
        }
    }

    Real CdsCmCouponPricer::swapletRate() const {
        return swapletPrice()/(accrualPeriod_ * discount_);
    }

    Rate CdsCmCouponPricer::capletRate(const Rate effectiveCap) const {
        return capletPrice(effectiveCap)/(accrualPeriod_ * discount_);
    }


    void CdsCmCouponPricer::initialize(const CmCdsCoupon& coupon) {
        Date today = Settings::instance().evaluationDate();
        coupon_ = &coupon;

        const Handle<DefaultProbabilityTermStructure>& defTS = 
            coupon_->creditIndex()->defaultProbTermStructure();
        const Handle<YieldTermStructure>& yTS =
            coupon_->creditIndex()->nominalTermStructure();

        // alternatively take the conventional recovery for the index
        // seniority; issue a warning rather?
        QL_REQUIRE( (coupon_->creditIndex()->seniority() == 
            recovery_->seniority()) || 
            recovery_->seniority() == NoSeniority,
            "Incompatible recovery quote seniority");
        fixingDate_    = coupon.fixingDate();
        refMaturity_   = fixingDate_ + coupon.creditIndex()->tenor();
        gearing_       = coupon.gearing();
        spread_        = coupon.spread();
        // no default event values:
        paymentDate_   = coupon.date();
        accrualPeriod_ = coupon.accrualPeriod();
        refProtStart_  = coupon.protectionStart();

        if(paymentDate_ > today) discount_ = 
            coupon.creditIndex()->nominalTermStructure()->discount(
                paymentDate_);
        else 
            discount_ = 1.;

		/* below is:
			effectiveProtectStart_ = 
			 max(today, refProtStart_, coupon.accrualStartDate())
		*/
        effectiveProtectStart_ = refProtStart_ <= today ? 
            today : refProtStart_;
		effectiveProtectStart_ = 
			coupon.accrualStartDate() <= effectiveProtectStart_ ? 
            effectiveProtectStart_ : coupon.accrualStartDate();

		if(!coupon_->hasOccurred()) {
            try{
				defaultProb_ = 
					defTS->defaultProbability(effectiveProtectStart_,
					  coupon_->protectionEnd());
            }catch(std::exception& e){
                std::string problemo(e.what());
            }
            survivalProb_ = 1. - defaultProb_;
            //butnot:defTS->survivalProbability(coupon_->referencePeriodEnd());
        }else{
            defaultProb_  = 0.;
            survivalProb_ = 1.;
        }
        // determine applicable rate
        if(today < fixingDate_ ) {
            // Although this is a call to index forecastFixing one can not call
            //   directly coupon_->creditIndex()->fixing; that would be at the 
            //   risk of using the conventional recovery rather than the quoted
            //   one. 
            boost::shared_ptr<CreditDefaultSwap> cdsSwap = 
                coupon_->creditIndex()->underlyingSwap(fixingDate_);
            cdsSwap->setPricingEngine(boost::shared_ptr<PricingEngine>(
                new MidPointCdsEngine(defTS, recovery_->value(), yTS)));
            effectiveRate_ = cdsSwap->fairSpread();
        }else{
            //'true' since most of the times I am calling this from a simulated
            //   scenario computation. Should it go into a constructor flag????
            effectiveRate_ = coupon_->creditIndex()->fixing(fixingDate_, true);
        }
    }
}
