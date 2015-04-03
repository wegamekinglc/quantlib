/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 Ferdinando Ametrano
 Copyright (C) 2007 François du Vignaud
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

/*! \file flatsmilesection.hpp
    \brief Flat SmileSection
*/

#ifndef quantlib_flat_smile_section_hpp
#define quantlib_flat_smile_section_hpp

#include <ql/termstructures/volatility/smilesection.hpp>

namespace QuantLib {

template <class T> class FlatSmileSection_t : public SmileSection_t<T> {
  public:
    FlatSmileSection_t(const Date &d, T vol, const DayCounter &dc,
                       const Date &referenceDate = Date(),
                       T atmLevel = Null<T>());
    FlatSmileSection_t(Time exerciseTime, T vol, const DayCounter &dc,
                       T atmLevel = Null<T>());
    //! \name SmileSection interface
    //@{
    T minStrike() const;
    T maxStrike() const;
    T atmLevel() const;
    //@}
  protected:
    T volatilityImpl(T) const;

  private:
    T vol_;
    T atmLevel_;
};

typedef FlatSmileSection_t<Real> FlatSmileSection;

template <class T> inline T FlatSmileSection_t<T>::minStrike() const {
    return T(QL_MIN_REAL);
}

template <class T> inline T FlatSmileSection_t<T>::maxStrike() const {
    return T(QL_MAX_REAL);
}

template <class T> inline T FlatSmileSection_t<T>::atmLevel() const {
    return this->atmLevel_;
}

template <class T> inline T FlatSmileSection_t<T>::volatilityImpl(T) const {
    return this->vol_;
}

// implementation

template <class T>
FlatSmileSection_t<T>::FlatSmileSection_t(const Date &d, T vol,
                                          const DayCounter &dc,
                                          const Date &referenceDate, T atmLevel)
    : SmileSection_t<T>(d, dc, referenceDate), vol_(vol), atmLevel_(atmLevel) {}

template <class T>
FlatSmileSection_t<T>::FlatSmileSection_t(Time exerciseTime, T vol,
                                          const DayCounter &dc, T atmLevel)
    : SmileSection_t<T>(exerciseTime, dc), vol_(vol), atmLevel_(atmLevel) {}
}

#endif
