/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2006, 2014 Ferdinando Ametrano
 Copyright (C) 2006 François du Vignaud
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

/*! \file capfloor.hpp
    \brief cap and floor class
*/

#ifndef quantlib_instruments_capfloor_hpp
#define quantlib_instruments_capfloor_hpp

#include <ql/instruments/capfloorbase.hpp>
#include <ql/pricingengines/capfloor/blackcapfloorengine.hpp>
#include <ql/math/solvers1d/newtonsafe.hpp>

namespace QuantLib {

namespace {

template <class T>
ImpliedVolHelperCapFloor<T>::ImpliedVolHelperCapFloor(
    const CapFloor_t<T> &cap,
    const Handle<YieldTermStructure_t<T> > &discountCurve, T targetValue,
    T displacement)
    : discountCurve_(discountCurve), targetValue_(targetValue) {

    // set an implausible value, so that calculation is forced
    // at first ImpliedVolHelperCapFloor::operator()(T x) call
    vol_ = boost::shared_ptr<SimpleQuote_t<T> >(new SimpleQuote_t<T>(-1));
    Handle<Quote_t<T> > h(vol_);
    engine_ = boost::shared_ptr<PricingEngine>(new BlackCapFloorEngine_t<T>(
        discountCurve_, h, Actual365Fixed(), displacement));
    cap.setupArguments(engine_->getArguments());

    results_ = dynamic_cast<const Instrument::results *>(engine_->getResults());
}

template <class T> T ImpliedVolHelperCapFloor<T>::operator()(T x) const {
    if (x != vol_->value()) {
        vol_->setValue(x);
        engine_->calculate();
    }
    return results_->value - targetValue_;
}

template <class T> T ImpliedVolHelperCapFloor<T>::derivative(T x) const {
    if (x != vol_->value()) {
        vol_->setValue(x);
        engine_->calculate();
    }
    std::map<std::string, boost::any>::const_iterator vega_ =
        results_->additionalResults.find("vega");
    QL_REQUIRE(vega_ != results_->additionalResults.end(), "vega not provided");
    return boost::any_cast<T>(vega_->second);
}
} // empty namespace

} // namespace QuantLib

#endif
