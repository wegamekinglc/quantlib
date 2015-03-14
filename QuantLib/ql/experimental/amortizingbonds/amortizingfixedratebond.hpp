/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Simon Ibbotson
 Copyright (C) 2015 Cheng Li

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

/*! \file amortizingfixedratebond.hpp
    \brief amortizing fixed-rate bond
*/

#ifndef quantlib_amortizing_fixed_rate_bond_hpp
#define quantlib_amortizing_fixed_rate_bond_hpp

#include <ql/instruments/bond.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/cashflows/cashflowvectors.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/time/schedule.hpp>

namespace QuantLib {

    class Schedule;

    using std::pow;
    //! amortizing fixed-rate bond
    template<class T = Real>
    class AmortizingFixedRateBond_t : public Bond_t<T> {
      public:
        AmortizingFixedRateBond_t(
                          Natural settlementDays,
                          const std::vector<T>& notionals,
                          const Schedule& schedule,
                          const std::vector<T>& coupons,
                          const DayCounter& accrualDayCounter,
                          BusinessDayConvention paymentConvention = Following,
                          const Date& issueDate = Date());
        /*! Automatically generates a set of equal coupons, with an
            amortizing bond.  The coupons are equal and the accrual
            daycount is only used for quoting/settlement purposes -
            not for calculating the coupons.
        */
        AmortizingFixedRateBond_t(
                          Natural settlementDays,
                          const Calendar& calendar,
                          T faceAmount,
                          const Date& startDate,
                          const Period& bondTenor,
                          const Frequency& sinkingFrequency,
                          const T coupon,
                          const DayCounter& accrualDayCounter,
                          BusinessDayConvention paymentConvention = Following,
                          const Date& issueDate = Date());
        Frequency frequency() const { return frequency_; }
        const DayCounter& dayCounter() const { return dayCounter_; }
      protected:
        Frequency frequency_;
        DayCounter dayCounter_;
    };

    typedef AmortizingFixedRateBond_t<Real> AmortizingFixedRateBond;

    template<class T>
    AmortizingFixedRateBond_t<T>::AmortizingFixedRateBond_t(
                                      Natural settlementDays,
                                      const std::vector<T>& notionals,
                                      const Schedule& schedule,
                                      const std::vector<T>& coupons,
                                      const DayCounter& accrualDayCounter,
                                      BusinessDayConvention paymentConvention,
                                      const Date& issueDate)
    : Bond_t<T>(settlementDays, schedule.calendar(), issueDate),
      frequency_(schedule.tenor().frequency()),
      dayCounter_(accrualDayCounter) {

        this->maturityDate_ = schedule.endDate();

        this->cashflows_ = FixedRateLeg_t<T>(schedule)
            .withNotionals(notionals)
            .withCouponRates(coupons, accrualDayCounter)
            .withPaymentAdjustment(paymentConvention);

        this->addRedemptionsToCashflows();

        QL_ENSURE(!this->cashflows().empty(), "bond with no cashflows!");
    }

    namespace  {

        std::pair<Integer,Integer> daysMinMax(const Period& p) {
            switch (p.units()) {
              case Days:
                return std::make_pair(p.length(), p.length());
              case Weeks:
                return std::make_pair(7*p.length(), 7*p.length());
              case Months:
                return std::make_pair(28*p.length(), 31*p.length());
              case Years:
                return std::make_pair(365*p.length(), 366*p.length());
              default:
                QL_FAIL("unknown time unit (" << Integer(p.units()) << ")");
            }
        }

        bool isSubPeriod(const Period& subPeriod,
                         const Period& superPeriod,
                         Integer& numSubPeriods) {

            std::pair<Integer, Integer> superDays(daysMinMax(superPeriod));
            std::pair<Integer, Integer> subDays(daysMinMax(subPeriod));

            //obtain the approximate time ratio
            Real minPeriodRatio =
                ((Real)superDays.first)/((Real)subDays.second);
            Real maxPeriodRatio =
                ((Real)superDays.second)/((Real)subDays.first);
            Integer lowRatio = static_cast<Integer>(std::floor(minPeriodRatio));
            Integer highRatio = static_cast<Integer>(std::ceil(maxPeriodRatio));

            try {
                for(Integer i=lowRatio; i <= highRatio; ++i) {
                    Period testPeriod = subPeriod * i;
                    if(testPeriod == superPeriod) {
                        numSubPeriods = i;
                        return true;
                    }
                }
            } catch(Error e) {
                return false;
            }

            return false;
        }

        Schedule sinkingSchedule(const Date& startDate,
                                 const Period& maturityTenor,
                                 const Frequency& sinkingFrequency,
                                 const Calendar& paymentCalendar) {
            Period freqPeriod(sinkingFrequency);
            Date maturityDate(startDate + maturityTenor);
            Schedule retVal(startDate, maturityDate, freqPeriod,
                            paymentCalendar, Unadjusted, Unadjusted,
                            DateGeneration::Backward, false);
            return retVal;
        }

        template<class T>
        std::vector<T> sinkingNotionals(const Period& maturityTenor,
                                           const Frequency& sinkingFrequency,
                                           T couponRate,
                                           T initialNotional) {
            Period freqPeriod(sinkingFrequency);
            Integer nPeriods;
            QL_REQUIRE(isSubPeriod(freqPeriod, maturityTenor, nPeriods),
                       "Bond frequency is incompatible with the maturity tenor");

            std::vector<T> notionals(nPeriods+1);
            notionals.front() = initialNotional;
            T coupon = couponRate / static_cast<Real>(sinkingFrequency);
            T compoundedInterest = 1.0;
            T totalValue = pow(1.0+coupon, nPeriods);
            for(Size i = 0; i < (Size)nPeriods-1; ++i) {
                compoundedInterest *= (1.0 + coupon);
                Real currentNotional = 0.0;
                if(coupon < 1.0e-12) {
                    currentNotional =
                       initialNotional*(1.0 - (i+1.0)/nPeriods);
                }
                else {
                    currentNotional =
                       initialNotional*(compoundedInterest - (compoundedInterest-1.0)/(1.0 - 1.0/totalValue));
                }
                notionals[i+1] = currentNotional;
            }
            notionals.back() = 0.0;
            return notionals;
        }

    }

    template<class T>
    AmortizingFixedRateBond_t<T>::AmortizingFixedRateBond_t(
                                      Natural settlementDays,
                                      const Calendar& calendar,
                                      T initialFaceAmount,
                                      const Date& startDate,
                                      const Period& bondTenor,
                                      const Frequency& sinkingFrequency,
                                      const T coupon,
                                      const DayCounter& accrualDayCounter,
                                      BusinessDayConvention paymentConvention,
                                      const Date& issueDate)
    : Bond_t<T>(settlementDays, calendar, issueDate),
      frequency_(sinkingFrequency),
      dayCounter_(accrualDayCounter) {

        QL_REQUIRE(bondTenor.length() > 0,
                   "bond tenor must be positive. "
                   << bondTenor << " is not allowed.");
        this->maturityDate_ = startDate + bondTenor;

        this->cashflows_ =
            FixedRateLeg_t<T>(sinkingSchedule(startDate, bondTenor,
                                         sinkingFrequency, calendar))
            .withNotionals(sinkingNotionals<T>(bondTenor,
                                            sinkingFrequency, coupon,
                                            initialFaceAmount))
            .withCouponRates(coupon, accrualDayCounter)
            .withPaymentAdjustment(paymentConvention);

        this->addRedemptionsToCashflows();
    }

}

#endif
