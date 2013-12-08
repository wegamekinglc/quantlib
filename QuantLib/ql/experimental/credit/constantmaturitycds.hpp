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


#ifndef quantlib_constant_maturity_cds_hpp
#define quantlib_constant_maturity_cds_hpp

#include <ql/instrument.hpp>
#include <ql/cashflow.hpp>
#include <ql/default.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/schedule.hpp>

#include <ql/experimental/credit/creditcmcoupon.hpp>

namespace QuantLib {

    class Claim;
    class SingleNameCreditIndex;

    class ConstantMaturityCDS : public Instrument {
    public:
        class arguments;
        class results;
        class engine;

        ConstantMaturityCDS(
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
            const Date& protectionStart = Date(),
            const boost::shared_ptr<Claim>& claim =
                boost::shared_ptr<Claim>()
            );

        //inspectors:
        Protection::Side side() const;
        Real notional() const;
        Rate gearing() const;
        bool settlesAccrual() const;
        bool paysAtDefaultTime() const;
        const Leg& coupons() const;
        const boost::shared_ptr<SingleNameCreditIndex>& 
                creditIndex() const;

        bool isExpired() const;
        void setupExpired() const;

        void setupArguments(PricingEngine::arguments* args) const;
        void fetchResults(const PricingEngine::results*) const;
        //calculated:
        Real couponLegNPV() const;
        Real defaultLegNPV() const;
        Real fairGearing() const;

    protected:
        //! Inception date:
        const Date protectionStart_;
        Leg couponLeg_;

        Protection::Side side_;
        Real notional_;
        boost::shared_ptr<SingleNameCreditIndex> creditIndex_;
        BusinessDayConvention paymentConvention_;
        DayCounter dayCount_;
        Real gearing_;
        bool settlesAccrual_;
        bool paysAtDefaultTime_;
        boost::shared_ptr<Claim> claim_;
        // pricing results:
        mutable Real fairGearing_;
        mutable Real couponLegNPV_;
        mutable Real defaultLegNPV_;
    };


    class ConstantMaturityCDS::arguments
        : public virtual PricingEngine::arguments {
      public:
        arguments();
        Protection::Side side;
        Real notional;
        Leg leg;
        Real gearing;
        bool settlesAccrual;
        bool paysAtDefaultTime;
        boost::shared_ptr<Claim> claim;
        Date protectionStart;
        void validate() const;
    };

    class ConstantMaturityCDS::results : public Instrument::results {
      public:
        Real couponLegNPV;
        Real defaultLegNPV;
        Real fairGearingFactor;

        void reset();
    };

    class ConstantMaturityCDS::engine
        : public GenericEngine<ConstantMaturityCDS::arguments,
                               ConstantMaturityCDS::results> {};

    // inlines
    inline Protection::Side ConstantMaturityCDS::side() const {
        return side_;
    }

    inline Real ConstantMaturityCDS::notional() const {
        return notional_;
    }

    inline Rate ConstantMaturityCDS::gearing() const {
        return gearing_;
    }

    inline bool ConstantMaturityCDS::settlesAccrual() const {
        return settlesAccrual_;
    }

    inline bool ConstantMaturityCDS::paysAtDefaultTime() const {
        return paysAtDefaultTime_;
    }

    inline const Leg& ConstantMaturityCDS::coupons() const {
        return couponLeg_;
    }

    inline const boost::shared_ptr<SingleNameCreditIndex>& 
        ConstantMaturityCDS::creditIndex() const {
        return creditIndex_;
    }

}

#endif
