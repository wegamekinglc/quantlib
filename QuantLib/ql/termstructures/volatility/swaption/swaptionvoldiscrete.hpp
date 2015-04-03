/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
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

/*! \file swaptionvoldiscrete.hpp
    \brief Discretized swaption volatility
*/

#ifndef quantlib_swaption_volatility_discrete_h
#define quantlib_swaption_volatility_discrete_h

#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/utilities/dataformatters.hpp>

namespace QuantLib {

template <class T>
class SwaptionVolatilityDiscrete_t : public LazyObject,
                                     public SwaptionVolatilityStructure_t<T> {
  public:
    SwaptionVolatilityDiscrete_t(const std::vector<Period> &optionTenors,
                                 const std::vector<Period> &swapTenors,
                                 Natural settlementDays, const Calendar &cal,
                                 BusinessDayConvention bdc,
                                 const DayCounter &dc);
    SwaptionVolatilityDiscrete_t(const std::vector<Period> &optionTenors,
                                 const std::vector<Period> &swapTenors,
                                 const Date &referenceDate, const Calendar &cal,
                                 BusinessDayConvention bdc,
                                 const DayCounter &dc);
    SwaptionVolatilityDiscrete_t(const std::vector<Date> &optionDates,
                                 const std::vector<Period> &swapTenors,
                                 const Date &referenceDate, const Calendar &cal,
                                 BusinessDayConvention bdc,
                                 const DayCounter &dc);
    const std::vector<Period> &optionTenors() const;
    const std::vector<Date> &optionDates() const;
    const std::vector<Time> &optionTimes() const;
    const std::vector<Period> &swapTenors() const;
    const std::vector<Time> &swapLengths() const;
    //@}
    //! \name Observer interface
    //@{
    void update();
    //@}
    //! \name LazyObject interface
    //@{
    void performCalculations() const;
    //@}
  protected:
    Size nOptionTenors_;
    std::vector<Period> optionTenors_;
    mutable std::vector<Date> optionDates_;
    mutable std::vector<Time> optionTimes_;
    mutable std::vector<Real> optionDatesAsReal_;
    mutable Interpolation optionInterpolator_;

    Size nSwapTenors_;
    std::vector<Period> swapTenors_;
    mutable std::vector<Time> swapLengths_;
    mutable Date evaluationDate_;

  private:
    void checkOptionTenors() const;
    void checkOptionDates() const;
    void checkSwapTenors() const;
    void initializeOptionDatesAndTimes() const;
    void initializeOptionTimes() const;
    void initializeSwapLengths() const;
};

typedef SwaptionVolatilityDiscrete_t<Real> SwaptionVolatilityDiscrete;

// inline

template <class T>
inline const std::vector<Period> &
SwaptionVolatilityDiscrete_t<T>::optionTenors() const {
    return optionTenors_;
}

template <class T>
inline const std::vector<Date> &
SwaptionVolatilityDiscrete_t<T>::optionDates() const {
    return optionDates_;
}

template <class T>
inline const std::vector<Time> &
SwaptionVolatilityDiscrete_t<T>::optionTimes() const {
    return optionTimes_;
}

template <class T>
inline const std::vector<Period> &
SwaptionVolatilityDiscrete_t<T>::swapTenors() const {
    return swapTenors_;
}

template <class T>
inline const std::vector<Time> &
SwaptionVolatilityDiscrete_t<T>::swapLengths() const {
    return swapLengths_;
}

// implementation

template <class T>
SwaptionVolatilityDiscrete_t<T>::SwaptionVolatilityDiscrete_t(
    const std::vector<Period> &optionTenors,
    const std::vector<Period> &swapTenors, Natural settlementDays,
    const Calendar &cal, BusinessDayConvention bdc, const DayCounter &dc)
    : SwaptionVolatilityStructure(settlementDays, cal, bdc, dc),
      nOptionTenors_(optionTenors.size()), optionTenors_(optionTenors),
      optionDates_(nOptionTenors_), optionTimes_(nOptionTenors_),
      optionDatesAsReal_(nOptionTenors_), nSwapTenors_(swapTenors.size()),
      swapTenors_(swapTenors), swapLengths_(nSwapTenors_) {

    checkOptionTenors();
    initializeOptionDatesAndTimes();

    checkSwapTenors();
    initializeSwapLengths();

    optionInterpolator_ = LinearInterpolation(
        optionTimes_.begin(), optionTimes_.end(), optionDatesAsReal_.begin());
    optionInterpolator_.update();
    optionInterpolator_.enableExtrapolation();

    registerWith(Settings::instance().evaluationDate());
    evaluationDate_ = Settings::instance().evaluationDate();
}

template <class T>
SwaptionVolatilityDiscrete_t<T>::SwaptionVolatilityDiscrete_t(
    const std::vector<Period> &optionTenors,
    const std::vector<Period> &swapTenors, const Date &referenceDate,
    const Calendar &cal, BusinessDayConvention bdc, const DayCounter &dc)
    : SwaptionVolatilityStructure(referenceDate, cal, bdc, dc),
      nOptionTenors_(optionTenors.size()), optionTenors_(optionTenors),
      optionDates_(nOptionTenors_), optionTimes_(nOptionTenors_),
      optionDatesAsReal_(nOptionTenors_), nSwapTenors_(swapTenors.size()),
      swapTenors_(swapTenors), swapLengths_(nSwapTenors_) {

    checkOptionTenors();
    initializeOptionDatesAndTimes();

    checkSwapTenors();
    initializeSwapLengths();

    optionInterpolator_ = LinearInterpolation(
        optionTimes_.begin(), optionTimes_.end(), optionDatesAsReal_.begin());
    optionInterpolator_.update();
    optionInterpolator_.enableExtrapolation();
}

template <class T>
SwaptionVolatilityDiscrete_t<T>::SwaptionVolatilityDiscrete_t(
    const std::vector<Date> &optionDates, const std::vector<Period> &swapTenors,
    const Date &referenceDate, const Calendar &cal, BusinessDayConvention bdc,
    const DayCounter &dc)
    : SwaptionVolatilityStructure(referenceDate, cal, bdc, dc),
      nOptionTenors_(optionDates.size()), optionTenors_(nOptionTenors_),
      optionDates_(optionDates), optionTimes_(nOptionTenors_),
      optionDatesAsReal_(nOptionTenors_), nSwapTenors_(swapTenors.size()),
      swapTenors_(swapTenors), swapLengths_(nSwapTenors_) {

    checkOptionDates();
    initializeOptionTimes();

    checkSwapTenors();
    initializeSwapLengths();

    optionInterpolator_ = LinearInterpolation(
        optionTimes_.begin(), optionTimes_.end(), optionDatesAsReal_.begin());
    optionInterpolator_.update();
    optionInterpolator_.enableExtrapolation();
}

template <class T>
void SwaptionVolatilityDiscrete_t<T>::checkOptionDates() const {
    QL_REQUIRE(optionDates_[0] > this->referenceDate(),
               "first option date ("
                   << optionDates_[0]
                   << ") must be greater than reference date ("
                   << this->referenceDate() << ")");
    for (Size i = 1; i < nOptionTenors_; ++i) {
        QL_REQUIRE(optionDates_[i] > optionDates_[i - 1],
                   "non increasing option dates: "
                       << io::ordinal(i) << " is " << optionDates_[i - 1]
                       << ", " << io::ordinal(i + 1) << " is "
                       << optionDates_[i]);
    }
}

template <class T>
void SwaptionVolatilityDiscrete_t<T>::checkOptionTenors() const {
    QL_REQUIRE(optionTenors_[0] > 0 * Days, "first option tenor is negative ("
                                                << optionTenors_[0] << ")");
    for (Size i = 1; i < nOptionTenors_; ++i)
        QL_REQUIRE(optionTenors_[i] > optionTenors_[i - 1],
                   "non increasing option tenor: "
                       << io::ordinal(i) << " is " << optionTenors_[i - 1]
                       << ", " << io::ordinal(i + 1) << " is "
                       << optionTenors_[i]);
}

template <class T>
void SwaptionVolatilityDiscrete_t<T>::checkSwapTenors() const {
    QL_REQUIRE(swapTenors_[0] > 0 * Days, "first swap tenor is negative ("
                                              << swapTenors_[0] << ")");
    for (Size i = 1; i < nSwapTenors_; ++i)
        QL_REQUIRE(swapTenors_[i] > swapTenors_[i - 1],
                   "non increasing swap tenor: "
                       << io::ordinal(i) << " is " << swapTenors_[i - 1] << ", "
                       << io::ordinal(i + 1) << " is " << swapTenors_[i]);
}

template <class T>
void SwaptionVolatilityDiscrete_t<T>::initializeOptionDatesAndTimes() const {
    for (Size i = 0; i < nOptionTenors_; ++i) {
        optionDates_[i] = this->optionDateFromTenor(optionTenors_[i]);
        optionDatesAsReal_[i] =
            static_cast<Real>(optionDates_[i].serialNumber());
    }
    initializeOptionTimes();
}

template <class T>
void SwaptionVolatilityDiscrete_t<T>::initializeOptionTimes() const {
    for (Size i = 0; i < nOptionTenors_; ++i)
        optionTimes_[i] = this->timeFromReference(optionDates_[i]);
}

template <class T>
void SwaptionVolatilityDiscrete_t<T>::initializeSwapLengths() const {
    for (Size i = 0; i < nSwapTenors_; ++i)
        swapLengths_[i] = this->swapLength(swapTenors_[i]);
}

template <class T>
void SwaptionVolatilityDiscrete_t<T>::performCalculations() const {
    // recalculate dates if necessary...
    if (this->moving_) {
        Date d = Settings::instance().evaluationDate();
        if (evaluationDate_ != d) {
            evaluationDate_ = d;
            initializeOptionDatesAndTimes();
            initializeSwapLengths();
            optionInterpolator_.update();
        }
    }
}

template <class T> void SwaptionVolatilityDiscrete_t<T>::update() {
    TermStructure::update();
    LazyObject::update();
}
}

#endif
