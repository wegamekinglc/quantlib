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

/*! \file smilesection.hpp
    \brief Smile section base class
*/

#ifndef quantlib_sabr_smile_section_hpp
#define quantlib_sabr_smile_section_hpp

#include <ql/termstructures/volatility/smilesection.hpp>
#include <ql/termstructures/volatility/sabr.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <vector>

namespace QuantLib {

template <class T> class SabrSmileSection_t : public SmileSection_t<T> {
  public:
    SabrSmileSection_t(Time timeToExpiry, T forward,
                       const std::vector<T> &sabrParameters);
    SabrSmileSection_t(const Date &d, T forward,
                       const std::vector<T> &sabrParameters,
                       const DayCounter &dc = Actual365Fixed());
    T minStrike() const { return 0.0; }
    T maxStrike() const { return QL_MAX_REAL; }
    T atmLevel() const { return forward_; }

  protected:
    T varianceImpl(T strike) const;
    T volatilityImpl(T strike) const;

  private:
    T alpha_, beta_, nu_, rho_, forward_;
};

typedef SabrSmileSection_t<Real> SabrSmileSection;

// implementation

template <class T>
SabrSmileSection_t<T>::SabrSmileSection_t(Time timeToExpiry, T forward,
                                          const std::vector<T> &sabrParams)
    : SmileSection_t<T>(timeToExpiry), forward_(forward) {

    alpha_ = sabrParams[0];
    beta_ = sabrParams[1];
    nu_ = sabrParams[2];
    rho_ = sabrParams[3];

    QL_REQUIRE(forward_ > 0.0, "at the money forward rate must be "
                               "positive: "
                                   << io::rate(forward_) << " not allowed");
    validateSabrParameters(alpha_, beta_, nu_, rho_);
}

template <class T>
SabrSmileSection_t<T>::SabrSmileSection_t(const Date &d, T forward,
                                          const std::vector<T> &sabrParams,
                                          const DayCounter &dc)
    : SmileSection_t<T>(d, dc), forward_(forward) {

    alpha_ = sabrParams[0];
    beta_ = sabrParams[1];
    nu_ = sabrParams[2];
    rho_ = sabrParams[3];

    QL_REQUIRE(forward_ > 0.0, "at the money forward rate must be "
                               "positive: "
                                   << io::rate(forward_) << " not allowed");
    validateSabrParameters(alpha_, beta_, nu_, rho_);
}

template <class T> T SabrSmileSection_t<T>::varianceImpl(T strike) const {
    strike = QLFCT::max(0.00001, strike);
    T vol = unsafeSabrVolatility(strike, forward_, this->exerciseTime(), alpha_,
                                 beta_, nu_, rho_);
    return vol * vol * this->exerciseTime();
}

template <class T> T SabrSmileSection_t<T>::volatilityImpl(T strike) const {
    strike = QLFCT::max(0.00001, strike);
    return unsafeSabrVolatility(strike, forward_, this->exerciseTime(), alpha_, beta_,
                                nu_, rho_);
}
}

#endif
