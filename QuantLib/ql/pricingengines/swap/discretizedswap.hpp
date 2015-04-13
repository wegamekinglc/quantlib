/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2004, 2007 StatPro Italia srl
 Copyright (C) 2015 Peter Caspers

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

/*! \file discretizedswap.hpp
    \brief Discretized swap class
*/

#ifndef quantlib_discretized_swap_hpp
#define quantlib_discretized_swap_hpp

#include <ql/instruments/vanillaswap.hpp>
#include <ql/discretizedasset.hpp>

namespace QuantLib {

template <class T> class DiscretizedSwap_t : public DiscretizedAsset_t<T> {
  public:
    DiscretizedSwap_t(const typename VanillaSwap_t<T>::arguments &,
                      const Date &referenceDate, const DayCounter &dayCounter);
    void reset(Size size);
    std::vector<Time> mandatoryTimes() const;

  protected:
    void preAdjustValuesImpl();
    void postAdjustValuesImpl();

  private:
    typename VanillaSwap_t<T>::arguments arguments_;
    std::vector<Time> fixedResetTimes_;
    std::vector<Time> fixedPayTimes_;
    std::vector<Time> floatingResetTimes_;
    std::vector<Time> floatingPayTimes_;
};

typedef DiscretizedSwap_t<Real> DiscretizedSwap;

// implementation

template <class T>
DiscretizedSwap_t<T>::DiscretizedSwap_t(const typename VanillaSwap_t<T>::arguments &args,
                                        const Date &referenceDate,
                                        const DayCounter &dayCounter)
    : arguments_(args) {

    fixedResetTimes_.resize(args.fixedResetDates.size());
    for (Size i = 0; i < fixedResetTimes_.size(); ++i)
        fixedResetTimes_[i] =
            dayCounter.yearFraction(referenceDate, args.fixedResetDates[i]);

    fixedPayTimes_.resize(args.fixedPayDates.size());
    for (Size i = 0; i < fixedPayTimes_.size(); ++i)
        fixedPayTimes_[i] =
            dayCounter.yearFraction(referenceDate, args.fixedPayDates[i]);

    floatingResetTimes_.resize(args.floatingResetDates.size());
    for (Size i = 0; i < floatingResetTimes_.size(); ++i)
        floatingResetTimes_[i] =
            dayCounter.yearFraction(referenceDate, args.floatingResetDates[i]);

    floatingPayTimes_.resize(args.floatingPayDates.size());
    for (Size i = 0; i < floatingPayTimes_.size(); ++i)
        floatingPayTimes_[i] =
            dayCounter.yearFraction(referenceDate, args.floatingPayDates[i]);
}

template <class T> void DiscretizedSwap_t<T>::reset(Size size) {
    this->values_ = Array_t<T>(size, 0.0);
    this->adjustValues();
}

template <class T>
std::vector<Time> DiscretizedSwap_t<T>::mandatoryTimes() const {
    std::vector<Time> times;
    for (Size i = 0; i < fixedResetTimes_.size(); i++) {
        Time t = fixedResetTimes_[i];
        if (t >= 0.0)
            times.push_back(t);
    }
    for (Size i = 0; i < fixedPayTimes_.size(); i++) {
        Time t = fixedPayTimes_[i];
        if (t >= 0.0)
            times.push_back(t);
    }
    for (Size i = 0; i < floatingResetTimes_.size(); i++) {
        Time t = floatingResetTimes_[i];
        if (t >= 0.0)
            times.push_back(t);
    }
    for (Size i = 0; i < floatingPayTimes_.size(); i++) {
        Time t = floatingPayTimes_[i];
        if (t >= 0.0)
            times.push_back(t);
    }
    return times;
}

template <class T> void DiscretizedSwap_t<T>::preAdjustValuesImpl() {
    // floating payments
    for (Size i = 0; i < floatingResetTimes_.size(); i++) {
        Time t = floatingResetTimes_[i];
        if (t >= 0.0 && this->isOnTime(t)) {
            DiscretizedDiscountBond_t<T> bond;
            bond.initialize(this->method(), floatingPayTimes_[i]);
            bond.rollback(this->time_);

            T nominal = arguments_.nominal;
            Time T0 = arguments_.floatingAccrualTimes[i];
            T spread = arguments_.floatingSpreads[i];
            T accruedSpread = nominal * T0 * spread;
            for (Size j = 0; j < this->values_.size(); j++) {
                T coupon = nominal * (1.0 - bond.values()[j]) +
                           accruedSpread * bond.values()[j];
                if (arguments_.type == VanillaSwap_t<T>::Payer)
                    this->values_[j] += coupon;
                else
                    this->values_[j] -= coupon;
            }
        }
    }
    // fixed payments
    for (Size i = 0; i < fixedResetTimes_.size(); i++) {
        Time t = fixedResetTimes_[i];
        if (t >= 0.0 && this->isOnTime(t)) {
            DiscretizedDiscountBond_t<T> bond;
            bond.initialize(this->method(), fixedPayTimes_[i]);
            bond.rollback(this->time_);

            T fixedCoupon = arguments_.fixedCoupons[i];
            for (Size j = 0; j < this->values_.size(); j++) {
                T coupon = fixedCoupon * bond.values()[j];
                if (arguments_.type == VanillaSwap_t<T>::Payer)
                    this->values_[j] -= coupon;
                else
                    this->values_[j] += coupon;
            }
        }
    }
}

template <class T> void DiscretizedSwap_t<T>::postAdjustValuesImpl() {
    // fixed coupons whose reset time is in the past won't be managed
    // in preAdjustValues()
    for (Size i = 0; i < fixedPayTimes_.size(); i++) {
        Time t = fixedPayTimes_[i];
        Time reset = fixedResetTimes_[i];
        if (t >= 0.0 && this->isOnTime(t) && reset < 0.0) {
            T fixedCoupon = arguments_.fixedCoupons[i];
            if (arguments_.type == VanillaSwap_t<T>::Payer)
                this->values_ -= fixedCoupon;
            else
                this->values_ += fixedCoupon;
        }
    }
    // the same applies to floating payments whose rate is already fixed
    for (Size i = 0; i < floatingPayTimes_.size(); i++) {
        Time t = floatingPayTimes_[i];
        Time reset = floatingResetTimes_[i];
        if (t >= 0.0 && this->isOnTime(t) && reset < 0.0) {
            T currentFloatingCoupon = arguments_.floatingCoupons[i];
            QL_REQUIRE(currentFloatingCoupon != Null<T>(),
                       "current floating coupon not given");
            if (arguments_.type == VanillaSwap_t<T>::Payer)
                this->values_ += currentFloatingCoupon;
            else
                this->values_ -= currentFloatingCoupon;
        }
    }
}

} // namespace QuantLib

#endif
