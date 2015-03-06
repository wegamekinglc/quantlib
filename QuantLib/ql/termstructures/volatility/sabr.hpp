/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
 Copyright (C) 2006 Mario Pucci
 Copyright (C) 2006 StatPro Italia srl
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

/*! \file sabr.hpp
    \brief SABR functions
*/

#ifndef quantlib_sabr_hpp
#define quantlib_sabr_hpp

#include <ql/types.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/math/comparison.hpp>
#include <ql/errors.hpp>

namespace QuantLib {

template <class T>
T unsafeSabrVolatility(T strike, T forward, Time expiryTime, T alpha, T beta,
                       T nu, T rho);

template <class T>
T sabrVolatility(T strike, T forward, Time expiryTime, T alpha, T beta, T nu,
                 T rho);

template <class T> void validateSabrParameters(T alpha, T beta, T nu, T rho);

// implementation

template <class T>
T unsafeSabrVolatility(T strike, T forward, Time expiryTime, T alpha, T beta,
                       T nu, T rho) {
    const T oneMinusBeta = 1.0 - beta;
    const T A = QLFCT::pow(forward * strike, oneMinusBeta);
    const T sqrtA = QLFCT::sqrt(A);
    T logM;
    if (!close(forward, strike))
        logM = QLFCT::log(forward / strike);
    else {
        const T epsilon = (forward - strike) / strike;
        logM = epsilon - .5 * epsilon * epsilon;
    }
    const T z = (nu / alpha) * sqrtA * logM;
    const T B = 1.0 - 2.0 * rho * z + z * z;
    const T C = oneMinusBeta * oneMinusBeta * logM * logM;
    const T tmp = (QLFCT::sqrt(B) + z - rho) / (1.0 - rho);
    const T xx = QLFCT::log(tmp);
    const T D = sqrtA * (1.0 + C / 24.0 + C * C / 1920.0);
    const T d =
        1.0 +
        expiryTime * (oneMinusBeta * oneMinusBeta * alpha * alpha / (24.0 * A) +
                      0.25 * rho * beta * nu * alpha / sqrtA +
                      (2.0 - 3.0 * rho * rho) * (nu * nu / 24.0));

    T multiplier;
    // computations become precise enough if the square of z worth
    // slightly more than the precision machine (hence the m)
    static const T m = 10;
    if (QLFCT::abs(z * z) > QL_EPSILON * m)
        multiplier = z / xx;
    else {
        multiplier =
            1.0 - 0.5 * rho * z - (3.0 * rho * rho - 2.0) * z * z / 12.0;
    }
    return (alpha / D) * multiplier * d;
}

template <class T> void validateSabrParameters(T alpha, T beta, T nu, T rho) {
    QL_REQUIRE(alpha > 0.0, "alpha must be positive: " << alpha
                                                       << " not allowed");
    QL_REQUIRE(beta >= 0.0 && beta <= 1.0,
               "beta must be in (0.0, 1.0): " << beta << " not allowed");
    QL_REQUIRE(nu >= 0.0, "nu must be non negative: " << nu << " not allowed");
    QL_REQUIRE(rho * rho < 1.0,
               "rho square must be less than one: " << rho << " not allowed");
}

template <class T>
T sabrVolatility(T strike, T forward, Time expiryTime, T alpha, T beta, T nu,
                 T rho) {
    QL_REQUIRE(strike > 0.0, "strike must be positive: " << io::rate(strike)
                                                         << " not allowed");
    QL_REQUIRE(forward > 0.0, "at the money forward rate must be "
                              "positive: "
                                  << io::rate(forward) << " not allowed");
    QL_REQUIRE(expiryTime >= 0.0, "expiry time must be non-negative: "
                                      << expiryTime << " not allowed");
    validateSabrParameters(alpha, beta, nu, rho);
    return unsafeSabrVolatility(strike, forward, expiryTime, alpha, beta, nu,
                                rho);
}
}

#endif
