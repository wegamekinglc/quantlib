/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
 Copyright (C) 2006 François du Vignaud
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

/*! \file interpolatedsmilesection.hpp
    \brief Interpolated smile section class
*/

#ifndef quantlib_interpolated_smile_section_hpp
#define quantlib_interpolated_smile_section_hpp

#include <ql/termstructure.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>

namespace QuantLib {

template <template <class> class Interpolator, class T>
class InterpolatedSmileSection_t : public SmileSection_t<T>, public LazyObject {
  public:
    InterpolatedSmileSection_t(
        Time expiryTime, const std::vector<T> &strikes,
        const std::vector<Handle<Quote_t<T> > > &stdDevHandles,
        const Handle<Quote_t<T> > &atmLevel,
        const Interpolator<T> &interpolator = Interpolator<T>(),
        const DayCounter &dc = Actual365Fixed());
    InterpolatedSmileSection_t(
        Time expiryTime, const std::vector<T> &strikes,
        const std::vector<T> &stdDevs, T atmLevel,
        const Interpolator<T> &interpolator = Interpolator<T>(),
        const DayCounter &dc = Actual365Fixed());

    InterpolatedSmileSection_t(
        const Date &d, const std::vector<T> &strikes,
        const std::vector<Handle<Quote_t<T> > > &stdDevHandles,
        const Handle<Quote_t<T> > &atmLevel,
        const DayCounter &dc = Actual365Fixed(),
        const Interpolator<T> &interpolator = Interpolator<T>(),
        const Date &referenceDate = Date());
    InterpolatedSmileSection_t(
        const Date &d, const std::vector<T> &strikes,
        const std::vector<T> &stdDevs, T atmLevel,
        const DayCounter &dc = Actual365Fixed(),
        const Interpolator<T> &interpolator = Interpolator<T>(),
        const Date &referenceDate = Date());
    void performCalculations() const;
    T varianceImpl(T strike) const;
    T volatilityImpl(T strike) const;
    T minStrike() const { return strikes_.front(); }
    T maxStrike() const { return strikes_.back(); }
    virtual T atmLevel() const { return atmLevel_->value(); }
    void update();

  private:
    Real exerciseTimeSquareRoot_;
    std::vector<T> strikes_;
    std::vector<Handle<Quote_t<T> > > stdDevHandles_;
    Handle<Quote_t<T> > atmLevel_;
    mutable std::vector<T> vols_;
    mutable Interpolation_t<T> interpolation_;
};

template <template <class> class Interpolator> struct InterpolatedSmileSection {
    typedef InterpolatedSmileSection_t<Interpolator, Real> Type;
};

template <template <class> class Interpolator, class T>
InterpolatedSmileSection_t<Interpolator, T>::InterpolatedSmileSection_t(
    Time timeToExpiry, const std::vector<T> &strikes,
    const std::vector<Handle<Quote_t<T> > > &stdDevHandles,
    const Handle<Quote_t<T> > &atmLevel, const Interpolator<T> &interpolator,
    const DayCounter &dc)
    : SmileSection_t<T>(timeToExpiry, dc),
      exerciseTimeSquareRoot_(std::sqrt(this->exerciseTime())), strikes_(strikes),
      stdDevHandles_(stdDevHandles), atmLevel_(atmLevel),
      vols_(stdDevHandles.size()) {
    for (Size i = 0; i < stdDevHandles_.size(); ++i)
        LazyObject::registerWith(stdDevHandles_[i]);
    LazyObject::registerWith(atmLevel_);
    // check strikes!!!!!!!!!!!!!!!!!!!!
    interpolation_ = interpolator.interpolate(strikes_.begin(), strikes_.end(),
                                              vols_.begin());
}

template <template <class> class Interpolator, class T>
InterpolatedSmileSection_t<Interpolator, T>::InterpolatedSmileSection_t(
    Time timeToExpiry, const std::vector<T> &strikes,
    const std::vector<T> &stdDevs, T atmLevel,
    const Interpolator<T> &interpolator, const DayCounter &dc)
    : SmileSection_t<T>(timeToExpiry, dc),
      exerciseTimeSquareRoot_(std::sqrt(this->exerciseTime())), strikes_(strikes),
      stdDevHandles_(stdDevs.size()), vols_(stdDevs.size()) {
    // fill dummy handles to allow generic handle-based
    // computations later on
    for (Size i = 0; i < stdDevs.size(); ++i)
        stdDevHandles_[i] = Handle<Quote_t<T> >(
            boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(stdDevs[i])));
    atmLevel_ = Handle<Quote_t<T> >(
        boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(atmLevel)));
    // check strikes!!!!!!!!!!!!!!!!!!!!
    interpolation_ = interpolator.interpolate(strikes_.begin(), strikes_.end(),
                                              vols_.begin());
}

template <template <class> class Interpolator, class T>
InterpolatedSmileSection_t<Interpolator, T>::InterpolatedSmileSection_t(
    const Date &d, const std::vector<T> &strikes,
    const std::vector<Handle<Quote_t<T> > > &stdDevHandles,
    const Handle<Quote_t<T> > &atmLevel, const DayCounter &dc,
    const Interpolator<T> &interpolator, const Date &referenceDate)
    : SmileSection(d, dc, referenceDate),
      exerciseTimeSquareRoot_(std::sqrt(this->exerciseTime())), strikes_(strikes),
      stdDevHandles_(stdDevHandles), atmLevel_(atmLevel),
      vols_(stdDevHandles.size()) {
    for (Size i = 0; i < stdDevHandles_.size(); ++i)
        LazyObject::registerWith(stdDevHandles_[i]);
    LazyObject::registerWith(atmLevel_);
    // check strikes!!!!!!!!!!!!!!!!!!!!
    interpolation_ = interpolator.interpolate(strikes_.begin(), strikes_.end(),
                                              vols_.begin());
}

template <template <class> class Interpolator, class T>
InterpolatedSmileSection_t<Interpolator, T>::InterpolatedSmileSection_t(
    const Date &d, const std::vector<T> &strikes, const std::vector<T> &stdDevs,
    T atmLevel, const DayCounter &dc, const Interpolator<T> &interpolator,
    const Date &referenceDate)
    : SmileSection_t<T>(d, dc, referenceDate),
      exerciseTimeSquareRoot_(std::sqrt(this->exerciseTime())), strikes_(strikes),
      stdDevHandles_(stdDevs.size()), vols_(stdDevs.size()) {
    // fill dummy handles to allow generic handle-based
    // computations later on
    for (Size i = 0; i < stdDevs.size(); ++i)
        stdDevHandles_[i] = Handle<Quote_t<T> >(
            boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(stdDevs[i])));
    atmLevel_ = Handle<Quote_t<T> >(
        boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(atmLevel)));
    // check strikes!!!!!!!!!!!!!!!!!!!!
    interpolation_ = interpolator.interpolate(strikes_.begin(), strikes_.end(),
                                              vols_.begin());
}

template <template <class> class Interpolator, class T>
inline void
InterpolatedSmileSection_t<Interpolator, T>::performCalculations() const {
    for (Size i = 0; i < stdDevHandles_.size(); ++i)
        vols_[i] = stdDevHandles_[i]->value() / exerciseTimeSquareRoot_;
    interpolation_.update();
}

#ifndef __DOXYGEN__
template <template <class> class Interpolator, class T>
T InterpolatedSmileSection_t<Interpolator, T>::varianceImpl(T strike) const {
    calculate();
    T v = interpolation_(strike, true);
    return v * v * this->exerciseTime();
}

template <template <class> class Interpolator, class T>
T InterpolatedSmileSection_t<Interpolator, T>::volatilityImpl(
    T strike) const {
    calculate();
    return interpolation_(strike, true);
}

template <template <class> class Interpolator, class T>
void InterpolatedSmileSection_t<Interpolator, T>::update() {
    LazyObject::update();
    SmileSection::update();
}
#endif

} // namespace QuantLib

#endif
