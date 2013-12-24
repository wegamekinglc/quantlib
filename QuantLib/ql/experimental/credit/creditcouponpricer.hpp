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

/*! \file creditcouponpricer.hpp
    \brief Pricing of a credit CM swap.
*/

#ifndef quantlib_credit_coupon_pricer_hpp
#define quantlib_credit_coupon_pricer_hpp

#include <ql/cashflow.hpp>
#include <ql/quote.hpp>
#include <ql/experimental/credit/creditcmcoupon.hpp>
#include <ql/experimental/credit/recoveryratequote.hpp>

namespace QuantLib {
    
    class YieldTermStructure;

    class CreditCouponPricer: public virtual Observer, 
                              public virtual Observable {
    public:
        virtual ~CreditCouponPricer() {}
        //! defaultable swaplet price
        virtual Real swapletPrice() const = 0;
        //! defaultable swaplet rate
        virtual Rate swapletRate() const  = 0;
        //! defaultable caplets
        virtual Real capletPrice(Rate effectiveCap) const = 0;
        virtual Rate capletRate(Rate effectiveCap) const = 0;
        //@}
        //! \name Observer interface
        //@{
        void update(){notifyObservers();}
        //@}
    protected:
        virtual void initialize(const CmCdsCoupon& coupon) = 0;
    };


    // \to do: add method to reset the vol handle (and readjust registration)
    //         add method to retrieve the vol handle
    //         add a smile
    /*! Black fwd cds cm coupon pricer. Volatility is swaption volatility for 
        the convexity adjustment presented here.

		Default events are priced including thus jump to default.
        The coupon does not knock out if defaults take place outside 
        the protection period. It does if they occur during protection. 
        The protection period does not necessarily coincides with the 
        accrual period, as in the case of a full leg. Thus the survival 
        probability on which the coupon payment is contingent is not the 
        survival from today to payment date but the one over the 
        protection period.
    */
    class CdsCmCouponPricer: public CreditCouponPricer {
    public:
        CdsCmCouponPricer(const Handle<Quote>& vol,
            const Handle<RecoveryRateQuote>& recovery
            ) 
         : vol_(vol),
           recovery_(recovery)
        {
           registerWith(vol_);
           registerWith(recovery_);
        }
        CdsCmCouponPricer(const Handle<Quote>& vol,
                          Real recovery
            ) 
         : vol_(vol),
           recovery_(boost::shared_ptr<RecoveryRateQuote>(
            new RecoveryRateQuote(recovery)))
        {
           registerWith(vol_);
        }
        CdsCmCouponPricer(Real vol, Real recovery);
        //! \name CreditCouponpricer Interface
        //@{
        Real swapletPrice() const;
        Rate swapletRate() const;
        Real capletPrice(Rate capValue) const;
        Rate capletRate(Rate capValue) const;
        void initialize(const CmCdsCoupon& coupon);
        //@}
    protected:
        // arguments and results:
        const CmCdsCoupon* coupon_;
        Date paymentDate_;
        Date fixingDate_;
        Date effectiveProtectStart_;
        Date refProtStart_;
		// default prob within its own period
        Real defaultProb_;
        // fwd conditional survival
        Real survivalProb_;
        Rate effectiveRate_;
        Real gearing_;
        Real spread_;
        //! discount value at the payment date
        mutable Real discount_;
        Date refMaturity_;
        mutable Time accrualPeriod_;
        Handle<Quote> vol_;
        Handle<RecoveryRateQuote> recovery_;
    };

    // \todo coupon price setters

    // \todo add forecast fixing coupon flag.

}

#endif
