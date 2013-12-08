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

#ifndef quantlib_credit_cm_coupon_hpp
#define quantlib_credit_cm_coupon_hpp

#include <ql/handle.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/time/schedule.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/index.hpp>

namespace QuantLib {

    class SingleNameCreditIndex;
    class YieldTermStructure;
    class CdsCmCouponPricer;

    /*! Default event contingent coupon. Defaults related to the protection \
    dates have potential impact on the coupon. It is for derived classes to \
    implement in which way (knockability, notional amortization, etc...)

    Notice theres no event type(s) or reference entity/basket at this level.
        @param protectionStartDate The start of the default applicability \
            period. Defaults affect in a certain way (to be specified by \
            derived classes) to the coupon payment.
        @param protectionEndDate The end of the default applicability period.
    */
    // Better name.... CreditCoupon? , RiskyCoupon? 
    class DefaultableCoupon : public Coupon {
    public:
        DefaultableCoupon(const Date& paymentDate,
                          Real nominal,
                          const Date& accrualStartDate,
                          const Date& accrualEndDate,
                          const Date& protectionStartDate,
                          const Date& protectionEndDate,
                          const Date& refPeriodStart = Date(),
                          const Date& refPeriodEnd = Date());
        const Date& protectionStart() const {return protectionStart_;}
        const Date& protectionEnd() const {return protectionEnd_;}
    protected:
        const Date protectionStart_;
        const Date protectionEnd_;
    };


    //! Constant maturity single name credit spread swap coupon
    /*! Coupon subject to an single name default event knockability and \
        the same single name index fixing. It always pays in arrears and \
        never fixes in arrears. It might settle the accrual.
        \warning This class does not perform any date adjustment,
                 i.e., the start and end date passed upon construction
                 should be already rolled to a business day.

        \to do: Incorrect naming? In the future a cmcds _index_ coupon cant be 
        named like this.
    */
    class CmCdsCoupon : public DefaultableCoupon, 
                        public Observer {
	public:
        /*! 
        @param gearing         Also, participation rate.
        @param fixingDays      Relative to referencePeriod/protection start
        @param startProtection Defaults on this date apply.
        @param endProtection   Defaults on this date do not apply.
        */
		CmCdsCoupon(const Date& paymentDate,
                    Real nominal,
                    const Date& startDate,
                    const Date& endDate,
                    Natural fixingDays,
                    const boost::shared_ptr<SingleNameCreditIndex>& index,
                    const Date& startProtection = Date(),
                    const Date& endProtection = Date(),
                    Real gearing = 1.0,
                    Spread spread = 0.0,
                    Rate cap = Null<Rate>(),
                    bool paysAccruedOnDefault = false,
                    bool settlesAtDefTime = false,
                    const Date& refPeriodStart = Date(),
                    const Date& refPeriodEnd = Date(),
                    const DayCounter& dayCounter = DayCounter()
                    );
        //! \name CashFlow interface
        //@{
        Real amount() const { return rate() * accrualPeriod() * nominal(); }
        //@}
        //! \name Coupon interface
        //@{
        //! Affected coupon rate, including gearing and spread
        Rate rate() const;
        DayCounter dayCounter() const { return dayCounter_; }
        //! Accrued amount due at the given date if no default took place.
        Real accruedAmount(const Date&) const;
        //@}
        //! \name Observer interface
        //@{
        void update();
        //@}
        Real price(const Handle<YieldTermStructure>& discountingCurve) const;
        //! \name Inspectors
        //@{
        //! credit index
        const boost::shared_ptr<SingleNameCreditIndex>& creditIndex() const {
            return index_;
        }
        const boost::shared_ptr<SingleNameCreditIndex>& index() const {
        ////const boost::shared_ptr<Index>& index() const {
            return index_;
        }
        //! fixing days
        Natural fixingDays() const { return fixingDays_; }
        //! fixing date
        virtual Date fixingDate() const { return fixingDate_; }
        //! If true the coupon pays accrual upon default events falling 
        //  within its reference period.
        bool settlesAccrual() const { return paysAccruedOnDefault_; }
        //! If accrual is paid upon default knock-out it is done at default 
        //  time, if false accrual is paid at coupon pay date.
        bool settlesAtDefTime() const { return settlesAtDefTime_; }
        //! index gearing, i.e. multiplicative coefficient for the index
        Real gearing() const { return gearing_; }
        //! spread paid over the fixing of the underlying credit index
        Spread spread() const { return spread_; }
        Rate cap() const {return cap_;}
        bool isCapped() const {return isCapped_;}
        //! fixing of the underlying index
        virtual Rate indexFixing() const;
        //! convexity adjustment
        virtual Rate convexityAdjustment() const;
        //! convexity-adjusted fixing
        virtual Rate adjustedFixing() const;
        // whether or not the coupon fixes in arrears
        // keep it for template compatibility?
        bool isInArrears() const { return false; }
        //@}
        //! \name Visitability
        //@{
        virtual void accept(AcyclicVisitor&);
        //@}

        // if theres to be a generic credit coupon these need to go to the base
        //   credit coupon class with the base CreditCouponpricer type
        void setPricer(const boost::shared_ptr<CdsCmCouponPricer>&);
        boost::shared_ptr<CdsCmCouponPricer> pricer() const {
            return pricer_;
        }
      private:
        //! convexity adjustment for the given index fixing
        Rate convexityAdjustmentImpl(Rate fixing) const;
        boost::shared_ptr<SingleNameCreditIndex> index_;
        DayCounter dayCounter_;
        Natural fixingDays_;
        Date fixingDate_;
        Real gearing_;
        Spread spread_;
        bool isCapped_;
        Rate cap_;
        boost::shared_ptr<CdsCmCouponPricer> pricer_;
        bool paysAccruedOnDefault_;
        bool settlesAtDefTime_;
	};

    // inlines
    inline Rate CmCdsCoupon::convexityAdjustment() const {
        return (gearing() == 0.0 ? 0.0 : adjustedFixing()-indexFixing());
    }

    //! CMCDS Leg helper. Builds a standard cm credit leg with protection 
    //    periods matching those of the accruals, and optionally a global
    //    protection starting at the first coupon.
    // \to do allow for other protection schedules (e.g. imm dates) through
    //    a rule. Achieve this through a second schedule?
    class CmCdsLeg {
    public:
        CmCdsLeg(const Schedule& schedule, 
                 const boost::shared_ptr<SingleNameCreditIndex>& index
                 );
        CmCdsLeg& withNotionals(Real notional);
        CmCdsLeg& withNotionals(const std::vector<Real>& notionals);
        CmCdsLeg& withPaymentDayCounter(const DayCounter&);
        CmCdsLeg& withPaymentAdjustment(BusinessDayConvention);
        CmCdsLeg& withFixingDays(Natural fixingDays);
        CmCdsLeg& withFixingDays(const std::vector<Natural>& fixingDays);
        CmCdsLeg& withGearings(Real gearing);
        CmCdsLeg& withGearings(const std::vector<Real>& gearings);
        CmCdsLeg& withSpreads(Spread spread);
        CmCdsLeg& withSpreads(const std::vector<Spread>& spreads);

        /* WIP PP */
        CmCdsLeg& withCaps(Rate cap);
        CmCdsLeg& withCaps(const std::vector<Rate>& caps);
       
        CmCdsLeg& withAccrualSettlement(bool flag = true);
        CmCdsLeg& withAccrualAtDefault(bool flag = true);
        /*! Establishes if protection start date (the date at which defaults 
           on this coupon are to be considered) is the reference start date,
           if false it will be the global schedule start date
        */
        CmCdsLeg& withGlobalProtectionDate(bool flag = true);

        operator Leg() const;
      private:
        Schedule schedule_;
        boost::shared_ptr<SingleNameCreditIndex> cdsIndex_;
        std::vector<Real> notionals_;
        DayCounter paymentDayCounter_;
        BusinessDayConvention paymentAdjustment_;
        std::vector<Natural> fixingDays_;
        std::vector<Real> gearings_;
        std::vector<Spread> spreads_;

        /* WIP PP*/
        std::vector<Rate> caps_;

        bool paysAccrual_;
        bool accrdAtDefault_;
        bool globalProtection_;
    };


}

#endif
