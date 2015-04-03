/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Mario Pucci
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

/*! \file spreadedsmilesection.hpp
    \brief Spreaded SmileSection class
*/

#ifndef quantlib_spreaded_smile_section_hpp
#define quantlib_spreaded_smile_section_hpp

#include <ql/termstructures/volatility/smilesection.hpp>
#include <ql/handle.hpp>
#include <ql/quote.hpp>

namespace QuantLib {

template <class T> class SpreadedSmileSection_t : public SmileSection_t<T> {
  public:
    SpreadedSmileSection_t(const boost::shared_ptr<SmileSection_t<T> > &,
                           const Handle<Quote_t<T> > &spread);
    //! \name SmileSection interface
    //@{
    T minStrike() const;
    T maxStrike() const;
    T atmLevel() const;
    const Date &exerciseDate() const;
    Time exerciseTime() const;
    const DayCounter &dayCounter() const;
    const Date &referenceDate() const;
    //@}
    //! \name LazyObject interface
    //@{
    void update() { this->notifyObservers(); }
    //@}
  protected:
    T volatilityImpl(T strike) const;

  private:
    const boost::shared_ptr<SmileSection_t<T> > underlyingSection_;
    const Handle<Quote_t<T> > spread_;
};

template <class T> inline T SpreadedSmileSection_t<T>::minStrike() const {
    return underlyingSection_->minStrike();
}

template <class T> inline T SpreadedSmileSection_t<T>::maxStrike() const {
    return underlyingSection_->maxStrike();
}

template <class T> inline T SpreadedSmileSection_t<T>::atmLevel() const {
    return underlyingSection_->atmLevel();
}

template <class T>
inline const Date &SpreadedSmileSection_t<T>::exerciseDate() const {
    return underlyingSection_->exerciseDate();
}

template <class T> inline Time SpreadedSmileSection_t<T>::exerciseTime() const {
    return underlyingSection_->exerciseTime();
}

template <class T>
inline const DayCounter &SpreadedSmileSection_t<T>::dayCounter() const {
    return underlyingSection_->dayCounter();
}

template <class T>
inline const Date &SpreadedSmileSection_t<T>::referenceDate() const {
    return underlyingSection_->referenceDate();
}

typedef SpreadedSmileSection_t<Real> SpreadedSmileSection;

// implementation

template <class T>
SpreadedSmileSection_t<T>::SpreadedSmileSection_t(
    const boost::shared_ptr<SmileSection_t<T> > &underlyingSection,
    const Handle<Quote_t<T> > &spread)
    : underlyingSection_(underlyingSection), spread_(spread) {
    this->registerWith(underlyingSection_);
    this->registerWith(spread_);
}

template <class T> T SpreadedSmileSection_t<T>::volatilityImpl(T k) const {
    return underlyingSection_->volatility(k) + spread_->value();
}
}

#endif
