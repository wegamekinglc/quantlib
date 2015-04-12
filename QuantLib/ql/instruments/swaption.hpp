/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2006 Cristina Duminuco
 Copyright (C) 2006 Marco Bianchetti
 Copyright (C) 2007 StatPro Italia srl
 Copyright (C) 2014 Ferdinando Ametrano
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

/*! \file swaption.hpp
    \brief Swaption class
*/

#ifndef quantlib_instruments_swaption_hpp
#define quantlib_instruments_swaption_hpp

#include <ql/instruments/swaptionbase.hpp>
#include <ql/option.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>

namespace QuantLib {

namespace {

template <class T> class ImpliedVolHelperSwaption {
  public:
    ImpliedVolHelperSwaption(
        const Swaption_t<T> &,
        const Handle<YieldTermStructure_t<T> > &discountCurve, T targetValue,
        T displacement);
    T operator()(T x) const;
    T derivative(T x) const;

  private:
    boost::shared_ptr<PricingEngine> engine_;
    Handle<YieldTermStructure_t<T> > discountCurve_;
    T targetValue_;
    boost::shared_ptr<SimpleQuote_t<T> > vol_;
    const typename Instrument_t<T>::results *results_;
};

template <class T>
ImpliedVolHelperSwaption<T>::ImpliedVolHelperSwaption(
    const Swaption_t<T> &swaption,
    const Handle<YieldTermStructure_t<T> > &discountCurve, T targetValue,
    T displacement)
    : discountCurve_(discountCurve), targetValue_(targetValue) {

    // set an implausible value, so that calculation is forced
    // at first ImpliedVolHelperSwaption::operator()(T x) call
    vol_ = boost::shared_ptr<SimpleQuote_t<T> >(new SimpleQuote_t<T>(-1.0));
    Handle<Quote_t<T> > h(vol_);
    engine_ = boost::shared_ptr<PricingEngine>(new BlackSwaptionEngine_t<T>(
        discountCurve_, h, Actual365Fixed(), displacement));
    swaption.setupArguments(engine_->getArguments());

    this->results_ =
        dynamic_cast<const Instrument::results *>(engine_->getResults());
}

template <class T> T ImpliedVolHelperSwaption<T>::operator()(T x) const {
    if (x != vol_->value()) {
        vol_->setValue(x);
        engine_->calculate();
    }
    return this->results_->value - targetValue_;
}

template <class T> T ImpliedVolHelperSwaption<T>::derivative(T x) const {
    if (x != vol_->value()) {
        vol_->setValue(x);
        engine_->calculate();
    }
    std::map<std::string, boost::any>::const_iterator vega_ =
        this->results_->additionalResults.find("vega");
    QL_REQUIRE(vega_ != results_->additionalResults.end(), "vega not provided");
    return boost::any_cast<T>(vega_->second);
}
}

template <class T>
Swaption_t<T>::Swaption_t(const boost::shared_ptr<VanillaSwap_t<T> > &swap,
                          const boost::shared_ptr<Exercise> &exercise,
                          Settlement::Type delivery)
    : Option_t<T>(boost::shared_ptr<Payoff>(), exercise), swap_(swap),
      settlementType_(delivery) {
    this->registerWith(swap_);
}

template <class T> bool Swaption_t<T>::isExpired() const {
    return detail::simple_event(this->exercise_->dates().back()).hasOccurred();
}

template <class T>
void Swaption_t<T>::setupArguments(PricingEngine::arguments *args) const {

    swap_->setupArguments(args);

    Swaption_t<T>::arguments *arguments =
        dynamic_cast<Swaption_t<T>::arguments *>(args);

    QL_REQUIRE(arguments != 0, "wrong argument type");

    arguments->swap = swap_;
    arguments->settlementType = settlementType_;
    arguments->exercise = this->exercise_;
}

template <class T> void Swaption_t<T>::arguments::validate() const {
    VanillaSwap_t<T>::arguments::validate();
    QL_REQUIRE(swap, "vanilla swap not set");
    QL_REQUIRE(this->exercise, "exercise not set");
}

template <class T>
T Swaption_t<T>::impliedVolatility(T targetValue,
                                   const Handle<YieldTermStructure_t<T> > &d,
                                   T guess, T accuracy, Natural maxEvaluations,
                                   T minVol, T maxVol, T displacement) const {
    // calculate();
    QL_REQUIRE(!isExpired(), "instrument expired");

    ImpliedVolHelperSwaption<T> f(*this, d, targetValue, displacement);
    // Brent solver;
    NewtonSafe_t<T> solver;
    solver.setMaxEvaluations(maxEvaluations);
    return solver.solve(f, accuracy, guess, minVol, maxVol);
}
} // namespace QuantLib

#endif
