/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Ferdinando Ametrano
 Copyright (C) 2006, 2007 StatPro Italia srl
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

/*! \file swaptionconstantvol.hpp
    \brief Constant swaption volatility
*/

#ifndef quantlib_swaption_constant_volatility_hpp
#define quantlib_swaption_constant_volatility_hpp

#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/time/period.hpp>
#include <ql/termstructures/volatility/flatsmilesection.hpp>
#include <ql/quotes/simplequote.hpp>

namespace QuantLib {

//! Constant swaption volatility, no time-strike dependence
template <class T>
class ConstantSwaptionVolatility_t : public SwaptionVolatilityStructure_t<T> {
  public:
    //! floating reference date, floating market data
    ConstantSwaptionVolatility_t(Natural settlementDays, const Calendar &cal,
                                 BusinessDayConvention bdc,
                                 const Handle<Quote_t<T> > &volatility,
                                 const DayCounter &dc);
    //! fixed reference date, floating market data
    ConstantSwaptionVolatility_t(const Date &referenceDate, const Calendar &cal,
                                 BusinessDayConvention bdc,
                                 const Handle<Quote_t<T> > &volatility,
                                 const DayCounter &dc);
    //! floating reference date, fixed market data
    ConstantSwaptionVolatility_t(Natural settlementDays, const Calendar &cal,
                                 BusinessDayConvention bdc, T volatility,
                                 const DayCounter &dc);
    //! fixed reference date, fixed market data
    ConstantSwaptionVolatility_t(const Date &referenceDate, const Calendar &cal,
                                 BusinessDayConvention bdc, T volatility,
                                 const DayCounter &dc);
    //! \name TermStructure interface
    //@{
    Date maxDate() const;
    //@}
    //! \name VolatilityTermStructure interface
    //@{
    T minStrike() const;
    T maxStrike() const;
    //@}
    //! \name SwaptionVolatilityStructure_t<T> interface
    //@{
    const Period &maxSwapTenor() const;
    //@}
  protected:
    boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(const Date &, const Period &) const;
    boost::shared_ptr<SmileSection_t<T> > smileSectionImpl(Time, Time) const;
    T volatilityImpl(const Date &, const Period &, T) const;
    T volatilityImpl(Time, Time, T) const;

  private:
    Handle<Quote_t<T> > volatility_;
    Period maxSwapTenor_;
};

// inline definitions

template <class T>
inline Date ConstantSwaptionVolatility_t<T>::maxDate() const {
    return Date::maxDate();
}

template <class T> inline T ConstantSwaptionVolatility_t<T>::minStrike() const {
    return QL_MIN_REAL;
}

template <class T> inline T ConstantSwaptionVolatility_t<T>::maxStrike() const {
    return QL_MAX_REAL;
}

template <class T>
inline const Period &ConstantSwaptionVolatility_t<T>::maxSwapTenor() const {
    return maxSwapTenor_;
}

typedef ConstantSwaptionVolatility_t<Real> ConstantSwaptionVolatility;

// implementation

// floating reference date, floating market data
template <class T>
ConstantSwaptionVolatility_t<T>::ConstantSwaptionVolatility_t(
    Natural settlementDays, const Calendar &cal, BusinessDayConvention bdc,
    const Handle<Quote_t<T> > &vol, const DayCounter &dc)
    : SwaptionVolatilityStructure_t<T>(settlementDays, cal, bdc, dc),
      volatility_(vol), maxSwapTenor_(100 * Years) {
    this->registerWith(volatility_);
}

// fixed reference date, floating market data
template <class T>
ConstantSwaptionVolatility_t<T>::ConstantSwaptionVolatility_t(
    const Date &referenceDate, const Calendar &cal, BusinessDayConvention bdc,
    const Handle<Quote_t<T> > &vol, const DayCounter &dc)
    : SwaptionVolatilityStructure_t<T>(referenceDate, cal, bdc, dc),
      volatility_(vol), maxSwapTenor_(100 * Years) {
    this->registerWith(volatility_);
}

// floating reference date, fixed market data
template <class T>
ConstantSwaptionVolatility_t<T>::ConstantSwaptionVolatility_t(
    Natural settlementDays, const Calendar &cal, BusinessDayConvention bdc,
    T vol, const DayCounter &dc)
    : SwaptionVolatilityStructure_t<T>(settlementDays, cal, bdc, dc),
      volatility_(boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(vol))),
      maxSwapTenor_(100 * Years) {}

// fixed reference date, fixed market data
template <class T>
ConstantSwaptionVolatility_t<T>::ConstantSwaptionVolatility_t(
    const Date &referenceDate, const Calendar &cal, BusinessDayConvention bdc,
    T vol, const DayCounter &dc)
    : SwaptionVolatilityStructure_t<T>(referenceDate, cal, bdc, dc),
      volatility_(boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(vol))),
      maxSwapTenor_(100 * Years) {}

template <class T>
boost::shared_ptr<SmileSection_t<T> >
ConstantSwaptionVolatility_t<T>::smileSectionImpl(const Date &d,
                                                  const Period &) const {
    T atmVol = volatility_->value();
    return boost::shared_ptr<SmileSection_t<T> >(new FlatSmileSection_t<T>(
        d, atmVol, this->dayCounter(), this->referenceDate()));
}

template <class T>
boost::shared_ptr<SmileSection_t<T> >
ConstantSwaptionVolatility_t<T>::smileSectionImpl(Time optionTime, Time) const {
    T atmVol = volatility_->value();
    return boost::shared_ptr<SmileSection_t<T> >(
        new FlatSmileSection_t<T>(optionTime, atmVol, this->dayCounter()));
}

template <class T>
T ConstantSwaptionVolatility_t<T>::volatilityImpl(const Date &, const Period &,
                                                  T) const {
    return volatility_->value();
}

template <class T>
T ConstantSwaptionVolatility_t<T>::volatilityImpl(Time, Time, T) const {
    return volatility_->value();
}
}

#endif
