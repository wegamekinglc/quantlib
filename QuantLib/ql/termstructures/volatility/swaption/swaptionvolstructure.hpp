/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl
 Copyright (C) 2006, 2008 Ferdinando Ametrano
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

/*! \file swaptionvolstructure.hpp
    \brief Swaption volatility structure
*/

#ifndef quantlib_swaption_volatility_structure_hpp
#define quantlib_swaption_volatility_structure_hpp

#include <ql/termstructures/voltermstructure.hpp>
#include <ql/math/rounding.hpp>

namespace QuantLib {

template <class T> class SmileSection_t;

//! %Swaption-volatility structure
/*! This abstract class defines the interface of concrete swaption
    volatility structures which will be derived from this one.
*/

template <class T>
class SwaptionVolatilityStructure_t : public VolatilityTermStructure_t<T> {
  public:
    /*! \name Constructors
        See the TermStructure documentation for issues regarding
        constructors.
    */
    //@{
    /*! \warning term structures initialized by means of this
                 constructor must manage their own reference date
                 by overriding the referenceDate() method.
        \deprecated
    */
    QL_DEPRECATED
    SwaptionVolatilityStructure_t(const Calendar &calendar,
                                  BusinessDayConvention bdc,
                                  const DayCounter &dc = DayCounter());

    /*! \warning term structures initialized by means of this
                 constructor must manage their own reference date
                 by overriding the referenceDate() method.
    */
    SwaptionVolatilityStructure_t(BusinessDayConvention bdc,
                                  const DayCounter &dc = DayCounter());
    //! initialize with a fixed reference date
    SwaptionVolatilityStructure_t(const Date &referenceDate,
                                  const Calendar &calendar,
                                  BusinessDayConvention bdc,
                                  const DayCounter &dc = DayCounter());
    //! calculate the reference date based on the global evaluation date
    SwaptionVolatilityStructure_t(Natural settlementDays, const Calendar &,
                                  BusinessDayConvention bdc,
                                  const DayCounter &dc = DayCounter());
    //@}
    virtual ~SwaptionVolatilityStructure_t() {}
    //! \name Volatility, variance and smile
    //@{
    //! returns the volatility for a given option tenor and swap tenor
    T volatility(const Period &optionTenor, const Period &swapTenor, T strike,
                 bool extrapolate = false) const;
    //! returns the volatility for a given option date and swap tenor
    T volatility(const Date &optionDate, const Period &swapTenor, T strike,
                 bool extrapolate = false) const;
    //! returns the volatility for a given option time and swap tenor
    T volatility(Time optionTime, const Period &swapTenor, T strike,
                 bool extrapolate = false) const;
    //! returns the volatility for a given option tenor and swap length
    T volatility(const Period &optionTenor, Time swapLength, T strike,
                 bool extrapolate = false) const;
    //! returns the volatility for a given option date and swap length
    T volatility(const Date &optionDate, Time swapLength, T strike,
                 bool extrapolate = false) const;
    //! returns the volatility for a given option time and swap length
    T volatility(Time optionTime, Time swapLength, T strike,
                 bool extrapolate = false) const;

    //! returns the Black variance for a given option tenor and swap tenor
    T blackVariance(const Period &optionTenor, const Period &swapTenor,
                    T strike, bool extrapolate = false) const;
    //! returns the Black variance for a given option date and swap tenor
    T blackVariance(const Date &optionDate, const Period &swapTenor, T strike,
                    bool extrapolate = false) const;
    //! returns the Black variance for a given option time and swap tenor
    T blackVariance(Time optionTime, const Period &swapTenor, T strike,
                    bool extrapolate = false) const;
    //! returns the Black variance for a given option tenor and swap length
    T blackVariance(const Period &optionTenor, Time swapLength, T strike,
                    bool extrapolate = false) const;
    //! returns the Black variance for a given option date and swap length
    T blackVariance(const Date &optionDate, Time swapLength, T strike,
                    bool extrapolate = false) const;
    //! returns the Black variance for a given option time and swap length
    T blackVariance(Time optionTime, Time swapLength, T strike,
                    bool extrapolate = false) const;

    //! returns the smile for a given option tenor and swap tenor
    boost::shared_ptr<SmileSection_t<T> >
    smileSection(const Period &optionTenor, const Period &swapTenor,
                 bool extr = false) const;
    //! returns the smile for a given option date and swap tenor
    boost::shared_ptr<SmileSection_t<T> > smileSection(const Date &optionDate,
                                                       const Period &swapTenor,
                                                       bool extr = false) const;
    //! returns the smile for a given option time and swap tenor
    boost::shared_ptr<SmileSection_t<T> > smileSection(Time optionTime,
                                                       const Period &swapTenor,
                                                       bool extr = false) const;
    //! returns the smile for a given option tenor and swap length
    boost::shared_ptr<SmileSection_t<T> >
    smileSection(const Period &optionTenor, Time swapLength,
                 bool extr = false) const;
    //! returns the smile for a given option date and swap length
    boost::shared_ptr<SmileSection_t<T> > smileSection(const Date &optionDate,
                                                       Time swapLength,
                                                       bool extr = false) const;
    //! returns the smile for a given option time and swap length
    boost::shared_ptr<SmileSection_t<T> >
    smileSection(Time optionTime, Time swapLength, bool extr = false) const;
    //@}
    //! \name Limits
    //@{
    //! the largest length for which the term structure can return vols
    virtual const Period &maxSwapTenor() const = 0;
    //! the largest swapLength for which the term structure can return vols
    Time maxSwapLength() const;
    //@}
    //! implements the conversion between swap tenor and swap (time) length
    Time swapLength(const Period &swapTenor) const;
    //! implements the conversion between swap dates and swap (time) length
    Time swapLength(const Date &start, const Date &end) const;

  protected:
    virtual boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(const Date &optionDate, const Period &swapTenor) const;
    virtual boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(Time optionTime, Time swapLength) const = 0;
    virtual T volatilityImpl(const Date &optionDate, const Period &swapTenor,
                             T strike) const;
    virtual T volatilityImpl(Time optionTime, Time swapLength,
                             T strike) const = 0;
    void checkSwapTenor(const Period &swapTenor, bool extrapolate) const;
    void checkSwapTenor(Time swapLength, bool extrapolate) const;
};

typedef SwaptionVolatilityStructure_t<Real> SwaptionVolatilityStructure;

// inline definitions

// 1. methods with Period-denominated exercise convert Period to Date and then
//    use the equivalent Date-denominated exercise methods
template <class T>
inline T
SwaptionVolatilityStructure_t<T>::volatility(const Period &optionTenor,
                                             const Period &swapTenor, T strike,
                                             bool extrapolate) const {
    Date optionDate = this->optionDateFromTenor(optionTenor);
    return volatility(optionDate, swapTenor, strike, extrapolate);
}

template <class T>
inline T
SwaptionVolatilityStructure_t<T>::volatility(const Period &optionTenor,
                                             Time swapLength, T strike,
                                             bool extrapolate) const {
    Date optionDate = this->optionDateFromTenor(optionTenor);
    return volatility(optionDate, swapLength, strike, extrapolate);
}

template <class T>
inline T SwaptionVolatilityStructure_t<T>::blackVariance(
    const Period &optionTenor, const Period &swapTenor, T strike,
    bool extrapolate) const {
    Date optionDate = this->optionDateFromTenor(optionTenor);
    return blackVariance(optionDate, swapTenor, strike, extrapolate);
}

template <class T>
inline T
SwaptionVolatilityStructure_t<T>::blackVariance(const Period &optionTenor,
                                                Time swapLength, T strike,
                                                bool extrapolate) const {
    Date optionDate = this->optionDateFromTenor(optionTenor);
    return blackVariance(optionDate, swapLength, strike, extrapolate);
}

template <class T>
inline boost::shared_ptr<SmileSection_t<T> >
SwaptionVolatilityStructure_t<T>::smileSection(const Period &optionTenor,
                                               const Period &swapTenor,
                                               bool extrapolate) const {
    Date optionDate = this->optionDateFromTenor(optionTenor);
    return smileSection(optionDate, swapTenor, extrapolate);
}

// 2. blackVariance methods rely on volatility methods
template <class T>
inline T SwaptionVolatilityStructure_t<T>::blackVariance(
    const Date &optionDate, const Period &swapTenor, T strike,
    bool extrapolate) const {
    T v = volatility(optionDate, swapTenor, strike, extrapolate);
    Time optionTime = this->timeFromReference(optionDate);
    return v * v * optionTime;
}

template <class T>
inline T SwaptionVolatilityStructure_t<T>::blackVariance(
    Time optionTime, const Period &swapTenor, T strike,
    bool extrapolate) const {
    T v = volatility(optionTime, swapTenor, strike, extrapolate);
    return v * v * optionTime;
}

template <class T>
inline T SwaptionVolatilityStructure_t<T>::blackVariance(
    const Date &optionDate, Time swapLength, T strike, bool extrapolate) const {
    T v = volatility(optionDate, swapLength, strike, extrapolate);
    Time optionTime = this->timeFromReference(optionDate);
    return v * v * optionTime;
}

template <class T>
inline T SwaptionVolatilityStructure_t<T>::blackVariance(
    Time optionTime, Time swapLength, T strike, bool extrapolate) const {
    T v = volatility(optionTime, swapLength, strike, extrapolate);
    return v * v * optionTime;
}

// 3. relying on xxxImpl methods
template <class T>
inline T
SwaptionVolatilityStructure_t<T>::volatility(const Date &optionDate,
                                             const Period &swapTenor, T strike,
                                             bool extrapolate) const {
    this->checkSwapTenor(swapTenor, extrapolate);
    this->checkRange(optionDate, extrapolate);
    this->checkStrike(strike, extrapolate);
    return volatilityImpl(optionDate, swapTenor, strike);
}

template <class T>
inline T SwaptionVolatilityStructure_t<T>::volatility(
    const Date &optionDate, Time swapLength, T strike, bool extrapolate) const {
    this->checkSwapTenor(swapLength, extrapolate);
    this->checkRange(optionDate, extrapolate);
    this->checkStrike(strike, extrapolate);
    Time optionTime = this->timeFromReference(optionDate);
    return volatilityImpl(optionTime, swapLength, strike);
}

template <class T>
inline T
SwaptionVolatilityStructure_t<T>::volatility(Time optionTime,
                                             const Period &swapTenor, T strike,
                                             bool extrapolate) const {
    this->checkSwapTenor(swapTenor, extrapolate);
    this->checkRange(optionTime, extrapolate);
    checkStrike(strike, extrapolate);
    Time length = swapLength(swapTenor);
    return volatilityImpl(optionTime, length, strike);
}

template <class T>
inline T
SwaptionVolatilityStructure_t<T>::volatility(Time optionTime, Time swapLength,
                                             T strike, bool extrapolate) const {
    this->checkSwapTenor(swapLength, extrapolate);
    this->checkRange(optionTime, extrapolate);
    this->checkStrike(strike, extrapolate);
    return volatilityImpl(optionTime, swapLength, strike);
}

template <class T>
inline boost::shared_ptr<SmileSection_t<T> >
SwaptionVolatilityStructure_t<T>::smileSection(const Date &optionDate,
                                               const Period &swapTenor,
                                               bool extrapolate) const {
    this->checkSwapTenor(swapTenor, extrapolate);
    this->checkRange(optionDate, extrapolate);
    return smileSectionImpl(optionDate, swapTenor);
}

template <class T>
inline boost::shared_ptr<SmileSection_t<T> >
SwaptionVolatilityStructure_t<T>::smileSection(Time optionTime, Time swapLength,
                                               bool extrapolate) const {
    this->checkSwapTenor(swapLength, extrapolate);
    this->checkRange(optionTime, extrapolate);
    return smileSectionImpl(optionTime, swapLength);
}

// 4. default implementation of Date-based xxxImpl methods
//    relying on the equivalent Time-based methods
template <class T>
inline boost::shared_ptr<SmileSection_t<T> >
SwaptionVolatilityStructure_t<T>::smileSectionImpl(const Date &optionDate,
                                                   const Period &swapT) const {
    return smileSectionImpl(this->timeFromReference(optionDate), swapLength(swapT));
}

template <class T>
inline T SwaptionVolatilityStructure_t<T>::volatilityImpl(
    const Date &optionDate, const Period &swapTenor, T strike) const {
    return volatilityImpl(this->timeFromReference(optionDate), swapLength(swapTenor),
                          strike);
}

template <class T>
inline Time SwaptionVolatilityStructure_t<T>::maxSwapLength() const {
    return swapLength(maxSwapTenor());
}

// implementation

template <class T>
SwaptionVolatilityStructure_t<T>::SwaptionVolatilityStructure_t(
    const Calendar &cal, BusinessDayConvention bdc, const DayCounter &dc)
    : VolatilityTermStructure_t<T>(bdc, dc) {
    this->calendar_ = cal;
}

template <class T>
SwaptionVolatilityStructure_t<T>::SwaptionVolatilityStructure_t(
    BusinessDayConvention bdc, const DayCounter &dc)
    : VolatilityTermStructure_t<T>(bdc, dc) {}

template <class T>
SwaptionVolatilityStructure_t<T>::SwaptionVolatilityStructure_t(
    const Date &referenceDate, const Calendar &calendar,
    BusinessDayConvention bdc, const DayCounter &dc)
    : VolatilityTermStructure_t<T>(referenceDate, calendar, bdc, dc) {}

template <class T>
SwaptionVolatilityStructure_t<T>::SwaptionVolatilityStructure_t(
    Natural settlementDays, const Calendar &calendar, BusinessDayConvention bdc,
    const DayCounter &dc)
    : VolatilityTermStructure_t<T>(settlementDays, calendar, bdc, dc) {}

template <class T>
Time SwaptionVolatilityStructure_t<T>::swapLength(const Period &p) const {
    QL_REQUIRE(p.length() > 0, "non-positive swap tenor (" << p << ") given");
    switch (p.units()) {
    case Months:
        return p.length() / 12.0;
    case Years:
        return static_cast<Time>(p.length());
    default:
        QL_FAIL("invalid Time Unit (" << p.units() << ") for swap length");
    }
}

template <class T>
Time SwaptionVolatilityStructure_t<T>::swapLength(const Date &start,
                                                  const Date &end) const {
    QL_REQUIRE(end > start, "swap end date ("
                                << end << ") must be greater than start ("
                                << start << ")");
    Time result = (end - start) / 365.25 * 12.0; // month unit
    result = ClosestRounding(0)(result);
    result /= 12.0; // year unit
    return result;
}

template <class T>
void SwaptionVolatilityStructure_t<T>::checkSwapTenor(const Period &swapTenor,
                                                      bool extrapolate) const {
    QL_REQUIRE(swapTenor.length() > 0, "non-positive swap tenor ("
                                           << swapTenor << ") given");
    QL_REQUIRE(extrapolate || this->allowsExtrapolation() ||
                   swapTenor <= maxSwapTenor(),
               "swap tenor (" << swapTenor << ") is past max tenor ("
                              << maxSwapTenor() << ")");
}

template <class T>
void SwaptionVolatilityStructure_t<T>::checkSwapTenor(Time swapLength,
                                                      bool extrapolate) const {
    QL_REQUIRE(swapLength > 0.0, "non-positive swap length (" << swapLength
                                                              << ") given");
    QL_REQUIRE(extrapolate || this->allowsExtrapolation() ||
                   swapLength <= maxSwapLength(),
               "swap tenor (" << swapLength << ") is past max tenor ("
                              << maxSwapLength() << ")");
}
}

#endif
