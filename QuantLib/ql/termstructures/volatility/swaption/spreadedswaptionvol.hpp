/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Ferdinando Ametrano
 Copyright (C) 2007 Giorgio Facchinetti
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

/*! \file spreadedswaptionvol.hpp
    \brief Spreaded swaption volatility
*/

#ifndef quantlib_spreaded_swaption_volstructure_h
#define quantlib_spreaded_swaption_volstructure_h

#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/termstructures/volatility/spreadedsmilesection.hpp>
#include <ql/quote.hpp>

namespace QuantLib {

template <class T>
class SpreadedSwaptionVolatility_t : public SwaptionVolatilityStructure {
  public:
    SpreadedSwaptionVolatility_t(const Handle<SwaptionVolatilityStructure> &,
                                 const Handle<Quote_t<T> > &spread);
    // All virtual methods of base classes must be forwarded
    //! \name TermStructure interface
    //@{
    DayCounter dayCounter() const;
    Date maxDate() const;
    Time maxTime() const;
    const Date &referenceDate() const;
    Calendar calendar() const;
    Natural settlementDays() const;
    //@}
    //! \name VolatilityTermStructure interface
    //@{
    T minStrike() const;
    T maxStrike() const;
    //@}
    //! \name SwaptionVolatilityStructure interface
    //@{
    const Period &maxSwapTenor() const;
    //@}
  protected:
    //! \name SwaptionVolatilityStructure interface
    //@{
    boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(const Date &optionDate, const Period &swapTenor) const;
    boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(Time optionTime, Time swapLength) const;
    T volatilityImpl(const Date &optionDate, const Period &swapTenor,
                     T strike) const;
    T volatilityImpl(Time optionTime, Time swapLength, T strike) const;
    //@}
  private:
    const Handle<SwaptionVolatilityStructure> baseVol_;
    const Handle<Quote_t<T> > spread_;
};

template <class T>
inline DayCounter SpreadedSwaptionVolatility_t<T>::dayCounter() const {
    return baseVol_->dayCounter();
}

template <class T>
inline Date SpreadedSwaptionVolatility_t<T>::maxDate() const {
    return baseVol_->maxDate();
}

template <class T>
inline Time SpreadedSwaptionVolatility_t<T>::maxTime() const {
    return baseVol_->maxTime();
}

template <class T>
inline const Date &SpreadedSwaptionVolatility_t<T>::referenceDate() const {
    return baseVol_->referenceDate();
}

template <class T>
inline Calendar SpreadedSwaptionVolatility_t<T>::calendar() const {
    return baseVol_->calendar();
}

template <class T>
inline Natural SpreadedSwaptionVolatility_t<T>::settlementDays() const {
    return baseVol_->settlementDays();
}

template <class T> inline T SpreadedSwaptionVolatility_t<T>::minStrike() const {
    return baseVol_->minStrike();
}

template <class T> inline T SpreadedSwaptionVolatility_t<T>::maxStrike() const {
    return baseVol_->maxStrike();
}

template <class T>
inline const Period &SpreadedSwaptionVolatility_t<T>::maxSwapTenor() const {
    return baseVol_->maxSwapTenor();
}

// implementation

template <class T>
SpreadedSwaptionVolatility_t<T>::SpreadedSwaptionVolatility_t(
    const Handle<SwaptionVolatilityStructure> &baseVol,
    const Handle<Quote_t<T> > &spread)
    : SwaptionVolatilityStructure(baseVol->businessDayConvention(),
                                  baseVol->dayCounter()),
      baseVol_(baseVol), spread_(spread) {
    enableExtrapolation(baseVol->allowsExtrapolation());
    registerWith(baseVol_);
    registerWith(spread_);
}

template <class T>
boost::shared_ptr<SmileSection_t<T> >
SpreadedSwaptionVolatility_t<T>::smileSectionImpl(const Date &d,
                                                  const Period &swapT) const {
    boost::shared_ptr<SmileSection_t<T> > baseSmile =
        baseVol_->smileSection(d, swapT, true);
    return boost::shared_ptr<SmileSection_t<T> >(
        new SpreadedSmileSection(baseSmile, spread_));
}

template <class T>
boost::shared_ptr<SmileSection_t<T> >
SpreadedSwaptionVolatility_t<T>::smileSectionImpl(Time optionTime,
                                                  Time swapLength) const {
    boost::shared_ptr<SmileSection_t<T> > baseSmile =
        baseVol_->smileSection(optionTime, swapLength, true);
    return boost::shared_ptr<SmileSection_t<T> >(
        new SpreadedSmileSection(baseSmile, spread_));
}

template <class T>
T SpreadedSwaptionVolatility_t<T>::volatilityImpl(const Date &d,
                                                  const Period &p,
                                                  T strike) const {
    return baseVol_->volatility(d, p, strike, true) + spread_->value();
}

template <class T>
T SpreadedSwaptionVolatility_t<T>::volatilityImpl(Time t, Time l,
                                                  T strike) const {
    return baseVol_->volatility(t, l, strike, true) + spread_->value();
}
}

#endif
