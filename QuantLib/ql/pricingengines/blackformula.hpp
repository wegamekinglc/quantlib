/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2003, 2004, 2005, 2006, 2008 Ferdinando Ametrano
 Copyright (C) 2006 Mark Joshi
 Copyright (C) 2006 StatPro Italia srl
 Copyright (C) 2007 Cristina Duminuco
 Copyright (C) 2007 Chiara Fornarola
 Copyright (C) 2013 Gary Kennedy
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

/*! \file blackformula.hpp
    \brief Black formula
*/

#ifndef quantlib_blackformula_hpp
#define quantlib_blackformula_hpp

#include <ql/option.hpp>
#include <ql/instruments/payoffs.hpp>
#include <ql/math/errorfunction.hpp>
#include <ql/math/solvers1d/newtonsafe.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif
#include <boost/math/special_functions/atanh.hpp>
#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic pop
#endif

namespace {
template <class T> void checkParameters(T strike, T forward, T displacement) {
    QL_REQUIRE(displacement >= T(0.0),
               "displacement (" << displacement << ") must be non-negative");
    QL_REQUIRE(strike + displacement >= T(0.0),
               "strike + displacement (" << strike << " + " << displacement
                                         << ") must be non-negative");
    QL_REQUIRE(forward + displacement > T(0.0),
               "forward + displacement (" << forward << " + " << displacement
                                          << ") must be positive");
}
}

namespace QuantLib {

/*! Black 1976 formula
    \warning instead of volatility it uses standard deviation,
             i.e. volatility*sqrt(timeToMaturity)
*/
template <class T>
T blackFormula(Option::Type optionType, T strike, T forward, T stdDev,
               T discount = 1.0, T displacement = 0.0) {
    checkParameters(strike, forward, displacement);
    QL_REQUIRE(stdDev >= T(0.0), "stdDev (" << stdDev
                                            << ") must be non-negative");
    QL_REQUIRE(discount > T(0.0), "discount (" << discount
                                               << ") must be positive");

    if (stdDev == T(0.0))
        return QLFCT::max<T>(T((forward - strike) * optionType), T(0.0)) *
               discount;
    forward = forward + displacement;
    strike = strike + displacement;

    // since displacement is non-negative strike==0 iff displacement==0
    // so returning forward*discount is OK
    if (strike == T(0.0))
        return (optionType == Option::Call ? forward * discount : T(0.0));

    T d1 = QLFCT::log(forward / strike) / stdDev + 0.5 * stdDev;
    T d2 = d1 - stdDev;
    CumulativeNormalDistribution_t<T> phi;
    T nd1 = phi(optionType * d1);
    T nd2 = phi(optionType * d2);
    T result = discount * optionType * (forward * nd1 - strike * nd2);
    QL_ENSURE(result >= T(0.0),
              "negative value (" << result << ") for " << stdDev << " stdDev, "
                                 << optionType << " option, " << strike
                                 << " strike , " << forward << " forward");
    return result;
}

/*! Black 1976 formula
    \warning instead of volatility it uses standard deviation,
             i.e. volatility*sqrt(timeToMaturity)
*/
template <class T>
T blackFormula(const boost::shared_ptr<PlainVanillaPayoff> &payoff, T forward,
               T stdDev, T discount = 1.0, T displacement = 0.0) {
    return blackFormula<T>(payoff->optionType(), payoff->strike(), forward,
                           stdDev, discount, displacement);
}

/*! Approximated Black 1976 implied standard deviation,
    i.e. volatility*sqrt(timeToMaturity).

    It is calculated using Brenner and Subrahmanyan (1988) and Feinstein
    (1988) approximation for at-the-money forward option, with the
    extended moneyness approximation by Corrado and Miller (1996)
    tape safe implementation
*/

template <class T>
T blackFormulaImpliedStdDevApproximation(Option::Type optionType, T strike,
                                         T forward, T blackPrice,
                                         T discount = 1.0,
                                         T displacement = 0.0) {
    checkParameters(strike, forward, displacement);
    QL_REQUIRE(blackPrice >= T(0.0),
               "blackPrice (" << blackPrice << ") must be non-negative");
    QL_REQUIRE(discount > T(0.0), "discount (" << discount
                                               << ") must be positive");

    forward = forward + displacement;
    strike = strike + displacement;

    // Brenner-Subrahmanyan (1988) and Feinstein (1988) ATM approx.
    T result0 = blackPrice / discount * QLFCT::sqrt(2.0 * M_PI) / forward;

    // Corrado and Miller extended moneyness approximation
    T moneynessDelta = optionType * (forward - strike);
    T moneynessDelta_2 = moneynessDelta / 2.0;
    T temp = blackPrice / discount - moneynessDelta_2;
    T moneynessDelta_PI = moneynessDelta * moneynessDelta / M_PI;
    T temp2 = temp * temp - moneynessDelta_PI;
    // approximation breaks down, 2 alternatives:
    // 1. zero it
    temp2 = QLFCT::CondExpLt(temp2, T(0.0), T(0.0), temp2);
    // 2. Manaster-Koehler (1982) efficient Newton-Raphson seed
    // return std::fabs(std::log(forward/strike))*std::sqrt(2.0);
    temp2 = QLFCT::sqrt(temp2);
    temp += temp2;
    temp *= QLFCT::sqrt(2.0 * M_PI);
    T result1 = temp / (forward + strike);

    T stdDev = QLFCT::CondExpEq(strike, forward, result0, result1);
    QL_ENSURE(stdDev >= T(0.0), "stdDev (" << stdDev
                                           << ") must be non-negative");
    return stdDev;
}

/*! Approximated Black 1976 implied standard deviation,
    i.e. volatility*sqrt(timeToMaturity).

    It is calculated using Brenner and Subrahmanyan (1988) and Feinstein
    (1988) approximation for at-the-money forward option, with the
    extended moneyness approximation by Corrado and Miller (1996)
    tape safe implementation
*/
template <class T>
T blackFormulaImpliedStdDevApproximation(
    const boost::shared_ptr<PlainVanillaPayoff> &payoff, T forward,
    T blackPrice, T discount = 1.0, T displacement = 0.0) {
    return blackFormulaImpliedStdDevApproximation<T>(
        payoff->optionType(), payoff->strike(), forward, blackPrice, discount,
        displacement);
}

// helper class
template <class T> class BlackImpliedStdDevHelper_t {
  public:
    BlackImpliedStdDevHelper_t(Option::Type optionType, T strike, T forward,
                               T undiscountedBlackPrice, T displacement = 0.0)
        : halfOptionType_(0.5 * optionType),
          signedStrike_(optionType * (strike + displacement)),
          signedForward_(optionType * (forward + displacement)),
          undiscountedBlackPrice_(undiscountedBlackPrice) {
        checkParameters(strike, forward, displacement);
        QL_REQUIRE(undiscountedBlackPrice >= T(0.0),
                   "undiscounted Black price (" << undiscountedBlackPrice
                                                << ") must be non-negative");
        signedMoneyness_ = optionType * QLFCT::log((forward + displacement) /
                                                   (strike + displacement));
    }
    T operator()(T stdDev) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(stdDev >= T(0.0), "stdDev (" << stdDev
                                                << ") must be non-negative");
#endif
        if (stdDev == T(0.0))
            return QLFCT::max<T>(signedForward_ - signedStrike_, T(0.0)) -
                   undiscountedBlackPrice_;
        T temp = halfOptionType_ * stdDev;
        T d = signedMoneyness_ / stdDev;
        T signedD1 = d + temp;
        T signedD2 = d - temp;
        T result = signedForward_ * N_(signedD1) - signedStrike_ * N_(signedD2);
        // numerical inaccuracies can yield a negative answer
        return QLFCT::max<T>(T(0.0), result) - undiscountedBlackPrice_;
    }
    T derivative(T stdDev) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(stdDev >= T(0.0), "stdDev (" << stdDev
                                                << ") must be non-negative");
#endif
        T signedD1 = signedMoneyness_ / stdDev + halfOptionType_ * stdDev;
        return signedForward_ * N_.derivative(signedD1);
    }

  private:
    T halfOptionType_;
    T signedStrike_, signedForward_;
    T undiscountedBlackPrice_, signedMoneyness_;
    CumulativeNormalDistribution_t<T> N_;
};

typedef BlackImpliedStdDevHelper_t<Real> BlackImpliedStdDevHelper;

/*! Black 1976 implied standard deviation,
    i.e. volatility*sqrt(timeToMaturity)
*/
template <class T>
T blackFormulaImpliedStdDev(Option::Type optionType, T strike, T forward,
                            T blackPrice, T discount = 1.0,
                            T displacement = 0.0, T guess = Null<T>(),
                            T accuracy = 1.0e-6, Natural maxIterations = 100) {
    checkParameters(strike, forward, displacement);

    QL_REQUIRE(discount > T(0.0), "discount (" << discount
                                               << ") must be positive");

    QL_REQUIRE(blackPrice >= T(0.0),
               "option price (" << blackPrice << ") must be non-negative");
    // check the price of the "other" option implied by put-call paity
    T otherOptionPrice =
        blackPrice - optionType * (forward - strike) * discount;
    QL_REQUIRE(
        otherOptionPrice >= T(0.0),
        "negative " << Option::Type(-1 * optionType) << " price ("
                    << otherOptionPrice
                    << ") implied by put-call parity. No solution exists for "
                    << optionType << " strike " << strike << ", forward "
                    << forward << ", price " << blackPrice << ", deflator "
                    << discount);

    // solve for the out-of-the-money option which has
    // greater vega/price ratio, i.e.
    // it is numerically more robust for implied vol calculations
    if (optionType == Option::Put && strike > forward) {
        optionType = Option::Call;
        blackPrice = otherOptionPrice;
    }
    if (optionType == Option::Call && strike < forward) {
        optionType = Option::Put;
        blackPrice = otherOptionPrice;
    }

    strike = strike + displacement;
    forward = forward + displacement;

    if (guess == Null<T>())
        guess = blackFormulaImpliedStdDevApproximation(
            optionType, strike, forward, blackPrice, discount, displacement);
    else
        QL_REQUIRE(guess >= T(0.0), "stdDev guess ("
                                        << guess << ") must be non-negative");
    BlackImpliedStdDevHelper f(optionType, strike, forward,
                               blackPrice / discount);
    NewtonSafe_t<T> solver;
    solver.setMaxEvaluations(maxIterations);
    T minSdtDev = T(0.0), maxStdDev = 24.0; // 24 = 300% * sqrt(60)
    T stdDev = solver.solve(f, accuracy, guess, minSdtDev, maxStdDev);
    QL_ENSURE(stdDev >= T(0.0), "stdDev (" << stdDev
                                           << ") must be non-negative");
    return stdDev;
}

/*! Black 1976 implied standard deviation,
    i.e. volatility*sqrt(timeToMaturity)
*/
template <class T>
T blackFormulaImpliedStdDev(const boost::shared_ptr<PlainVanillaPayoff> &payoff,
                            T forward, T blackPrice, T discount = 1.0,
                            T displacement = 0.0, T guess = Null<T>(),
                            T accuracy = 1.0e-6, Natural maxIterations = 100) {
    return blackFormulaImpliedStdDev(
        payoff->optionType(), payoff->strike(), forward, blackPrice, discount,
        displacement, guess, accuracy, maxIterations);
}

/*! Black 1976 probability of being in the money (in the bond martingale
    measure), i.e. N(d2).
    It is a risk-neutral probability, not the real world one.
    \warning instead of volatility it uses standard deviation,
             i.e. volatility*sqrt(timeToMaturity)
*/
template <class T>
T blackFormulaCashItmProbability(Option::Type optionType, T strike, T forward,
                                 T stdDev, T displacement = 0.0) {
    checkParameters(strike, forward, displacement);
    if (stdDev == T(0.0))
        return (forward * optionType > strike * optionType ? T(1.0) : T(0.0));

    forward = forward + displacement;
    strike = strike + displacement;
    if (strike == T(0.0))
        return (optionType == Option::Call ? T(1.0) : T(0.0));
    T d2 = QLFCT::log(forward / strike) / stdDev - 0.5 * stdDev;
    CumulativeNormalDistribution_t<T> phi;
    return phi(optionType * d2);
}

/*! Black 1976 probability of being in the money (in the bond martingale
    measure), i.e. N(d2).
    It is a risk-neutral probability, not the real world one.
    \warning instead of volatility it uses standard deviation,
             i.e. volatility*sqrt(timeToMaturity)
*/
template <class T>
T blackFormulaCashItmProbability(
    const boost::shared_ptr<PlainVanillaPayoff> &payoff, T forward, T stdDev,
    T displacement = 0.0) {
    return blackFormulaCashItmProbability(
        payoff->optionType(), payoff->strike(), forward, stdDev, displacement);
}

/*! Black 1976 formula for standard deviation derivative
    \warning instead of volatility it uses standard deviation, i.e.
             volatility*sqrt(timeToMaturity), and it returns the
             derivative with respect to the standard deviation.
             If T is the time to maturity Black vega would be
             blackStdDevDerivative(strike, forward, stdDev)*sqrt(T)
*/
template <class T>
T blackFormulaStdDevDerivative(T strike, T forward, T stdDev, T discount = 1.0,
                               T displacement = 0.0) {
    checkParameters(strike, forward, displacement);
    QL_REQUIRE(stdDev >= T(0.0), "stdDev (" << stdDev
                                            << ") must be non-negative");
    QL_REQUIRE(discount > T(0.0), "discount (" << discount
                                               << ") must be positive");

    forward = forward + displacement;
    strike = strike + displacement;

    if (stdDev == T(0.0) || strike == T(0.0))
        return T(0.0);

    T d1 = QLFCT::log(forward / strike) / stdDev + .5 * stdDev;
    return discount * forward *
           CumulativeNormalDistribution_t<T>().derivative(d1);
}

/*! Black 1976 formula for  derivative with respect to implied vol, this
    is basically the vega, but if you want 1% change multiply by 1%
*/
template <class T>
T blackFormulaVolDerivative(T strike, T forward, T stdDev, Time expiry,
                            T discount = 1.0, T displacement = 0.0) {
    return blackFormulaStdDevDerivative(strike, forward, stdDev, discount,
                                        displacement) *
           QLFCT::sqrt(expiry);
}

/*! Black 1976 formula for standard deviation derivative
    \warning instead of volatility it uses standard deviation, i.e.
             volatility*sqrt(timeToMaturity), and it returns the
             derivative with respect to the standard deviation.
             If T is the time to maturity Black vega would be
             blackStdDevDerivative(strike, forward, stdDev)*sqrt(T)
*/
template <class T>
T blackFormulaStdDevDerivative(
    const boost::shared_ptr<PlainVanillaPayoff> &payoff, T forward, T stdDev,
    T discount = 1.0, T displacement = 0.0) {
    return blackFormulaStdDevDerivative(payoff->strike(), forward, stdDev,
                                        discount, displacement);
}

/*! Black style formula when forward is normal rather than
    log-normal. This is essentially the model of Bachelier.

    \warning Bachelier model needs absolute volatility, not
             percentage volatility. Standard deviation is
             absoluteVolatility*sqrt(timeToMaturity)
*/
template <class T>
T bachelierBlackFormula(Option::Type optionType, T strike, T forward, T stdDev,
                        T discount = 1.0) {
    QL_REQUIRE(stdDev >= T(0.0), "stdDev (" << stdDev
                                            << ") must be non-negative");
    QL_REQUIRE(discount > T(0.0), "discount (" << discount
                                               << ") must be positive");
    T d = (forward - strike) * optionType, h = d / stdDev;
    if (stdDev == T(0.0))
        return discount * QLFCT::max(d, T(0.0));
    CumulativeNormalDistribution_t<T> phi;
    T result = discount * (stdDev * phi.derivative(h) + d * phi(h));
    QL_ENSURE(result >= T(0.0),
              "negative value (" << result << ") for " << stdDev << " stdDev, "
                                 << optionType << " option, " << strike
                                 << " strike , " << forward << " forward");
    return result;
}

template <class T> inline T h(T eta) {

    const static T A0 = 3.994961687345134e-1;
    const static T A1 = 2.100960795068497e+1;
    const static T A2 = 4.980340217855084e+1;
    const static T A3 = 5.988761102690991e+2;
    const static T A4 = 1.848489695437094e+3;
    const static T A5 = 6.106322407867059e+3;
    const static T A6 = 2.493415285349361e+4;
    const static T A7 = 1.266458051348246e+4;

    const static T B0 = 1.000000000000000e+0;
    const static T B1 = 4.990534153589422e+1;
    const static T B2 = 3.093573936743112e+1;
    const static T B3 = 1.495105008310999e+3;
    const static T B4 = 1.323614537899738e+3;
    const static T B5 = 1.598919697679745e+4;
    const static T B6 = 2.392008891720782e+4;
    const static T B7 = 3.608817108375034e+3;
    const static T B8 = -2.067719486400926e+2;
    const static T B9 = 1.174240599306013e+1;

    QL_REQUIRE(eta >= T(0.0), "eta (" << eta << ") must be non-negative");

    const T num =
        A0 +
        eta * (A1 +
               eta * (A2 +
                      eta * (A3 +
                             eta * (A4 + eta * (A5 + eta * (A6 + eta * A7))))));

    const T den =
        B0 +
        eta *
            (B1 +
             eta *
                 (B2 +
                  eta *
                      (B3 +
                       eta * (B4 +
                              eta * (B5 +
                                     eta * (B6 +
                                            eta * (B7 +
                                                   eta * (B8 + eta * B9))))))));

    return QLFCT::sqrt(eta) * (num / den);
}

/*! Black style formula when forward is normal rather than
    log-normal. This is essentially the model of Bachelier.

    \warning Bachelier model needs absolute volatility, not
             percentage volatility. Standard deviation is
             absoluteVolatility*sqrt(timeToMaturity)
*/
template <class T>
T bachelierBlackFormula(const boost::shared_ptr<PlainVanillaPayoff> &payoff,
                        T forward, T stdDev, T discount = 1.0) {
    return bachelierBlackFormula(payoff->optionType(), payoff->strike(),
                                 forward, stdDev, discount);
}

/*! Approximated Bachelier implied volatility

    It is calculated using  the analytic implied volatility approximation
    of J. Choi, K Kim and M. Kwak (2009), “Numerical Approximation of the
    Implied Volatility Under Arithmetic Brownian Motion”,
    Applied Math. Finance, 16(3), pp. 261-268.

    WARNING: This does not work yet with CppAD due to the missing atanh ...
*/
template <class T>
T bachelierBlackFormulaImpliedVol(Option::Type optionType, T strike, T forward,
                                  Time tte, T bachelierPrice, T discount = 1.0) {
    const static T SQRT_QL_EPSILON = QLFCT::sqrt(QL_EPSILON);

    QL_REQUIRE(tte > T(0.0), "tte (" << tte << ") must be positive");

    T forwardPremium = bachelierPrice / discount;

    T straddlePremium;
    if (optionType == Option::Call) {
        straddlePremium = 2.0 * forwardPremium - (forward - strike);
    } else {
        straddlePremium = 2.0 * forwardPremium + (forward - strike);
    }

    T nu = (forward - strike) / straddlePremium;
    QL_REQUIRE(nu <= T(1.0), "nu (" << nu << ") must be <= 1.0");
    QL_REQUIRE(nu >= T(-1.0), "nu (" << nu << ") must be >= -1.0");

    nu = QLFCT::max(T(-1.0 + QL_EPSILON), QLFCT::min(nu, T(1.0 - QL_EPSILON)));

    // nu / arctanh(nu) -> 1 as nu -> 0
    T eta =
        (QLFCT::abs(nu) < SQRT_QL_EPSILON) ? 1.0 : nu / boost::math::atanh(nu);

    T heta = h(eta);

    T impliedBpvol = QLFCT::sqrt(M_PI / (2 * tte)) * straddlePremium * heta;

    return impliedBpvol;
}
}

#endif
