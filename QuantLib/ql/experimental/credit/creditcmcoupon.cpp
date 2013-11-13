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

#include <ql/experimental/credit/creditcmcoupon.hpp>
#include <ql/cashflows/cashflowvectors.hpp>
#include <ql/experimental/credit/creditindex.hpp>
#include <ql/experimental/credit/creditcouponpricer.hpp>

namespace QuantLib {

    DefaultableCoupon::DefaultableCoupon(
        const Date &paymentDate, 
        Real nominal, 
        const Date& accrualStartDate, 
        const Date& accrualEndDate, 
        const Date& protectionStartDate, 
        const Date& protectionEndDate, 
        const Date& refPeriodStart, 
        const Date& refPeriodEnd)
    : Coupon(paymentDate, nominal, accrualStartDate, accrualEndDate, 
             refPeriodStart, refPeriodEnd),
      protectionStart_(protectionStartDate == Date() ? 
            accrualStartDate : protectionStartDate), 
      protectionEnd_(protectionEndDate == Date() ? 
            accrualEndDate : protectionEndDate) 
    {
        // Date constraints
        QL_REQUIRE(protectionStartDate < protectionEndDate, 
            "Incoherent default observation period dates in coupon");
        QL_REQUIRE(accrualStartDate_ < accrualEndDate_, 
            "Accrual dates seem in reverse order.");
        // equals because last protection date is not included in the p.period
        QL_REQUIRE(paymentDate_ >= protectionEnd_, 
            "Payment date lies before end of protection.");
    }

    void CmCdsCoupon::update() {
        if(pricer_) pricer_->initialize(*this);
        notifyObservers(); 
    }

    CmCdsCoupon::CmCdsCoupon(
        const Date& paymentDate,
        Real nominal,
        const Date& accrualStartDate, 
        const Date& accrualEndDate, 
        Natural fixingDays,
        const boost::shared_ptr<SingleNameCreditIndex>& 
           index,
        const Date& startProtection,
        const Date& endProtection,
        Real gearing,
        Spread spread,
        Rate cap,
        bool paysAccruedOnDefault,
        bool settlesAtDefTime,
        const Date& refPeriodStart, 
        const Date& refPeriodEnd, 
        const DayCounter& dayCounter) 
    : DefaultableCoupon(paymentDate, nominal, accrualStartDate, accrualEndDate, 
        startProtection, endProtection,
        refPeriodStart, refPeriodEnd),
      index_(index),
      dayCounter_(dayCounter),
      fixingDays_(fixingDays),
      fixingDate_(index->fixingCalendar().advance(accrualStartDate_,
          -static_cast<Integer>(fixingDays), Days, Preceding)),
      gearing_(gearing),
      spread_(spread),
      cap_(cap),    
      isCapped_(cap != Null<Rate>() ? true : false ),
      paysAccruedOnDefault_(paysAccruedOnDefault),
      settlesAtDefTime_(settlesAtDefTime)
    {
        // strict zero comparison because protection includes the start 
        //   date, this has to be coherent with the flags when calling 
        //   issuer.defaultedBetween(...)
        QL_REQUIRE(static_cast<Integer>(fixingDays) > 0, 
            "CMCDS coupon can not fix within the reference period");//not in arrears
        QL_REQUIRE(gearing_ > 0., "CMCDS gearing must be positive");


        registerWith(index_);
        registerWith(Settings::instance().evaluationDate());
    }

    // move to credit coupon base class..................
    Rate CmCdsCoupon::rate() const {
        QL_REQUIRE(pricer_, "pricer not set");
        pricer_->initialize(*this);
        return pricer_->swapletRate();// even with the cap now ????
    }

    void CmCdsCoupon::accept(AcyclicVisitor& v) {
        Visitor<CmCdsCoupon>* vcmc =
            dynamic_cast<Visitor<CmCdsCoupon>*>(&v);
        if (vcmc != 0)
            vcmc->visit(*this);
        else
            Coupon::accept(v);
    }

    // or inline it and add the extra dependency...
    Rate CmCdsCoupon::indexFixing() const {
        return index_->fixing(fixingDate());
    }

    Real CmCdsCoupon::price(
        const Handle<YieldTermStructure>& discountingCurve) const {
        return amount() * discountingCurve->discount(date());
    }

    Rate CmCdsCoupon::adjustedFixing() const {
        return (rate()-spread())/gearing();
    }

    void CmCdsCoupon::setPricer(
        const boost::shared_ptr<CdsCmCouponPricer>& pricer) {
        if (pricer_)
            unregisterWith(pricer_);
        pricer_ = pricer;
        if (pricer_)
            registerWith(pricer_);
        update();
    }

    Real CmCdsCoupon::accruedAmount(const Date& d) const {
        if (d <= accrualStartDate_ || d > paymentDate_) {
            return 0.0;
        } else {
            return nominal() * rate() * 
                dayCounter().yearFraction(accrualStartDate_,
                                          std::min(d, accrualEndDate_),
                                          refPeriodStart_,
                                          refPeriodEnd_);
        }
    }



    // CMCDS Leg helper definitions ------------------------------------------

    CmCdsLeg::CmCdsLeg(const Schedule& schedule,
                       const boost::shared_ptr<SingleNameCreditIndex>& index)
    : schedule_(schedule), cdsIndex_(index),
      paymentAdjustment_(Following),
      paysAccrual_(false),
      accrdAtDefault_(false),
      globalProtection_(false) {}

    CmCdsLeg& CmCdsLeg::withAccrualSettlement(bool flag) {
        paysAccrual_ = flag;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withAccrualAtDefault(bool flag) {
        accrdAtDefault_ = flag;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withNotionals(Real notional) {
        notionals_ = std::vector<Real>(1, notional);
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withNotionals(const std::vector<Real>& notionals) {
        notionals_ = notionals;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withPaymentDayCounter(const DayCounter& dayCounter) {
        paymentDayCounter_ = dayCounter;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withPaymentAdjustment(
        BusinessDayConvention convention) {
        paymentAdjustment_ = convention;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withFixingDays(Natural fixingDays) {
        fixingDays_ = std::vector<Natural>(1, fixingDays);
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withFixingDays(
        const std::vector<Natural>& fixingDays) {
        fixingDays_ = fixingDays;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withGearings(Real gearing) {
        gearings_ = std::vector<Real>(1, gearing);
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withGearings(const std::vector<Real>& gearings) {
        gearings_ = gearings;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withSpreads(Spread spread) {
        spreads_ = std::vector<Spread>(1, spread);
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withSpreads(const std::vector<Spread>& spreads) {
        spreads_ = spreads;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withCaps(Rate cap) {
        caps_ = std::vector<Rate>(1, cap);
        return *this;
    }
  
    CmCdsLeg& CmCdsLeg::withCaps(const std::vector<Rate>& caps) {
        caps_ = caps;
        return *this;
    }

    CmCdsLeg& CmCdsLeg::withGlobalProtectionDate(bool flag) {
        globalProtection_ = flag;
        return *this;
    }

    CmCdsLeg::operator Leg() const {
        Size n = schedule_.size()-1;
        QL_REQUIRE(!notionals_.empty(), "no notional given");
        QL_REQUIRE(notionals_.size() <= n,
                   "too many nominals (" << notionals_.size() <<
                   "), only " << n << " required");
        QL_REQUIRE(gearings_.size()<=n,
                   "too many gearings (" << gearings_.size() <<
                   "), only " << n << " required");
        QL_REQUIRE(spreads_.size()<=n,
                   "too many spreads (" << spreads_.size() <<
                   "), only " << n << " required");
        QL_REQUIRE(caps_.size()<=n,
                   "too many caps (" << caps_.size() <<
                   "), only " << n << " required");
	
        Leg leg; leg.reserve(n);

        Calendar calendar = schedule_.calendar();

        Date refStart, start, refEnd, end, protStart, protEnd;
        Date lastPaymentDate = calendar.adjust(schedule_.date(n), 
            paymentAdjustment_);

        for (Size i=0; i<n; ++i) {
            start = schedule_.date(i);
            if(globalProtection_) 
                protStart = schedule_.date(0);
            else 
                protStart = start;
            refStart = start;
            refEnd   = end = schedule_.date(i+1);
            Date paymentDate = calendar.adjust(end, paymentAdjustment_);
            if (i==n-1 && !schedule_.isRegular(i+1)) {
                BusinessDayConvention bdc = schedule_.businessDayConvention();
                refEnd = schedule_.calendar().adjust( 
                    start + schedule_.tenor(), bdc);
            }

            // Few things missing??? (caps..) check...  I am assuming here that the coupons are alll floaters

            protEnd = refEnd;
                leg.push_back(boost::shared_ptr<CashFlow>(new CmCdsCoupon(
                        paymentDate,
                        detail::get(notionals_, i, 1.0),
                        start, end,
                        detail::get(fixingDays_, i, 0),
                        cdsIndex_,
                        protStart, protEnd,
                        detail::get(gearings_, i, 1.0),
                        detail::get(spreads_, i, 0.0),
                        detail::get(caps_, i, Null<Rate>()),
                        paysAccrual_,
                        accrdAtDefault_,
                        refStart,
                        refEnd,
                        paymentDayCounter_)));
        }
        return leg;
    }

}
