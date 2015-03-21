/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2006 Cristina Duminuco
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

#ifndef quantlib_abcd_hpp
#define quantlib_abcd_hpp

#include <ql/types.hpp>
#include <ql/errors.hpp>
#include <ql/math/comparison.hpp>

namespace QuantLib {

template <class Ty>
inline void validateAbcdParameters(Ty a,
                                   Ty, // no condition on b
                                   Ty c, Ty d) {
    QL_REQUIRE(a + d >= 0, "a (" << a << ") + d (" << d
                                 << ") must be non negative");
    QL_REQUIRE(c >= 0, "c (" << c << ") must be non negative");
    QL_REQUIRE(d >= 0, "d (" << d << ") must be non negative");
}

//! %Abcd functional form for instantaneous volatility
/*! \f[ f(T-t) = [ a + b(T-t) ] e^{-c(T-t)} + d \f]
    following Rebonato's notation. */
template <class Ty> class AbcdFunction_t : public std::unary_function<Ty, Ty> {

  public:
    AbcdFunction_t(Ty a = -0.06, Ty b = 0.17, Ty c = 0.54, Ty d = 0.17);

    //! volatility function value at time u: \f[ f(u) \f]
    Ty operator()(Time u) const;

    //! time at which the volatility function reaches maximum (if any)
    Ty maximumLocation() const;

    //! maximum value of the volatility function
    Ty maximumVolatility() const;

    //! volatility function value at time 0: \f[ f(0) \f]
    Ty shortTermVolatility() const { return a_ + d_; }

    //! volatility function value at time +inf: \f[ f(\inf) \f]
    Ty longTermVolatility() const { return d_; }

    /*! instantaneous covariance function at time t between T-fixing and
        S-fixing rates \f[ f(T-t)f(S-t) \f] */
    Ty covariance(Time t, Time T, Time S) const;

    /*! integral of the instantaneous covariance function between
        time t1 and t2 for T-fixing and S-fixing rates
        \f[ \int_{t1}^{t2} f(T-t)f(S-t)dt \f] */
    Ty covariance(Time t1, Time t2, Time T, Time S) const;

    /*! average volatility in [tMin,tMax] of T-fixing rate:
       \f[ \sqrt{ \frac{\int_{tMin}^{tMax} f^2(T-u)du}{tMax-tMin} } \f] */
    Ty volatility(Time tMin, Time tMax, Time T) const;

    /*! variance between tMin and tMax of T-fixing rate:
        \f[ \frac{\int_{tMin}^{tMax} f^2(T-u)du}{tMax-tMin} \f] */
    Ty variance(Time tMin, Time tMax, Time T) const;

    // INSTANTANEOUS
    /*! instantaneous volatility at time t of the T-fixing rate:
        \f[ f(T-t) \f] */
    Ty instantaneousVolatility(Time t, Time T) const;

    /*! instantaneous variance at time t of T-fixing rate:
        \f[ f(T-t)f(T-t) \f] */
    Ty instantaneousVariance(Time t, Time T) const;

    /*! instantaneous covariance at time t between T and S fixing rates:
        \f[ f(T-u)f(S-u) \f] */
    Ty instantaneousCovariance(Time u, Time T, Time S) const;

    // PRIMITIVE
    /*! indefinite integral of the instantaneous covariance function at
        time t between T-fixing and S-fixing rates
        \f[ \int f(T-t)f(S-t)dt \f] */
    Ty primitive(Time t, Time T, Time S) const;

    /*! Inspectors */
    Ty a() const { return a_; }
    Ty b() const { return b_; }
    Ty c() const { return c_; }
    Ty d() const { return d_; }

  private:
    Ty a_, b_, c_, d_;
};

typedef AbcdFunction_t<Real> AbcdFunction;

// Helper class used by unit tests
template <class Ty> class AbcdSquared_t : public std::unary_function<Ty, Ty> {

  public:
    AbcdSquared_t(Ty a, Ty b, Ty c, Ty d, Time T, Time S);
    Ty operator()(Time t) const;

  private:
    boost::shared_ptr<AbcdFunction> abcd_;
    Time T_, S_;
};

typedef AbcdSquared_t<Real> AbcdSquared;

template <class Ty> inline Ty abcdBlackVolatility(Time u, Ty a, Ty b, Ty c, Ty d) {
    AbcdFunction_t<Ty> model(a, b, c, d);
    return model.volatility(0., u, u);
}

// implementation

template <class Ty>
AbcdFunction_t<Ty>::AbcdFunction_t(Ty a, Ty b, Ty c, Ty d)
    : a_(a), b_(b), c_(c), d_(d) {
    validateAbcdParameters(a, b, c, d);
}

template <class Ty> Ty AbcdFunction_t<Ty>::operator()(Time u) const {
    return u < 0 ? 0.0 : (a_ + b_ * u) * QLFCT::exp(-c_ * u) + d_;
}

template <class Ty> Ty AbcdFunction_t<Ty>::maximumLocation() const {
    if (b_ <= 0) {
        return 0.0;
    } else {
        if ((b_ - c_ * a_) / (c_ * b_) > 0) {
            return (b_ - c_ * a_) / (c_ * b_);
        } else
            return 0.0;
    }
}

template <class Ty> Ty AbcdFunction_t<Ty>::maximumVolatility() const {
    if (b_ <= 0) {
        return shortTermVolatility();
    } else {
        if ((b_ - c_ * a_) / (c_ * b_) > 0.0) {
            return b_ / c_ * QLFCT::exp(-1.0 + c_ * a_ / b_) + d_;
        } else
            return shortTermVolatility();
    }
}

template <class Ty>
Ty AbcdFunction_t<Ty>::volatility(Time tMin, Time tMax, Time T) const {
    if (tMax == tMin)
        return instantaneousVolatility(tMax, T);
    QL_REQUIRE(tMax > tMin, "tMax must be > tMin");
    return QLFCT::sqrt(variance(tMin, tMax, T) / (tMax - tMin));
}

template <class Ty>
Ty AbcdFunction_t<Ty>::variance(Time tMin, Time tMax, Time T) const {
    return covariance(tMin, tMax, T, T);
}

template <class Ty>
Ty AbcdFunction_t<Ty>::covariance(Time t, Time T, Time S) const {
    return (*this)(T - t) * (*this)(S - t);
}

template <class Ty>
Ty AbcdFunction_t<Ty>::covariance(Time t1, Time t2, Time T, Time S) const {
    QL_REQUIRE(t1 <= t2, "integrations bounds (" << t1 << "," << t2
                                                 << ") are in reverse order");
    Time cutOff = QLFCT::min(S, T);
    if (t1 >= cutOff) {
        return 0.0;
    } else {
        cutOff = QLFCT::min(t2, cutOff);
        return primitive(cutOff, T, S) - primitive(t1, T, S);
    }
}

// INSTANTANEOUS
template <class Ty>
Ty AbcdFunction_t<Ty>::instantaneousVolatility(Time u, Time T) const {
    return QLFCT::sqrt(instantaneousVariance(u, T));
}

template <class Ty>
Ty AbcdFunction_t<Ty>::instantaneousVariance(Time u, Time T) const {
    return instantaneousCovariance(u, T, T);
}
template <class Ty>
Ty AbcdFunction_t<Ty>::instantaneousCovariance(Time u, Time T, Time S) const {
    return (*this)(T - u) * (*this)(S - u);
}

// PRIMITIVE
template <class Ty>
Ty AbcdFunction_t<Ty>::primitive(Time t, Time T, Time S) const {
    if (T < t || S < t)
        return 0.0;

    if (close(c_, 0.0)) {
        Ty v = a_ + d_;
        return t *
               (v * v + v * b_ * S + v * b_ * T - v * b_ * t + b_ * b_ * S * T -
                0.5 * b_ * b_ * t * (S + T) + b_ * b_ * t * t / 3.0);
    }

    Ty k1 = QLFCT::exp(c_ * t), k2 = QLFCT::exp(c_ * S), k3 = QLFCT::exp(c_ * T);

    return (b_ * b_ * (-1 - 2 * c_ * c_ * S * T - c_ * (S + T) +
                       k1 * k1 * (1 + c_ * (S + T - 2 * t) +
                                  2 * c_ * c_ * (S - t) * (T - t))) +
            2 * c_ * c_ *
                (2 * d_ * a_ * (k2 + k3) * (k1 - 1) + a_ * a_ * (k1 * k1 - 1) +
                 2 * c_ * d_ * d_ * k2 * k3 * t) +
            2 * b_ * c_ * (a_ * (-1 - c_ * (S + T) +
                                 k1 * k1 * (1 + c_ * (S + T - 2 * t))) -
                           2 * d_ * (k3 * (1 + c_ * S) + k2 * (1 + c_ * T) -
                                     k1 * k3 * (1 + c_ * (S - t)) -
                                     k1 * k2 * (1 + c_ * (T - t))))) /
           (4 * c_ * c_ * c_ * k2 * k3);
}

//===========================================================================//
//                               AbcdSquared                                //
//===========================================================================//

template <class Ty>
AbcdSquared_t<Ty>::AbcdSquared_t(Ty a, Ty b, Ty c, Ty d, Time T, Time S)
    : abcd_(new AbcdFunction_t<Ty>(a, b, c, d)), T_(T), S_(S) {}

template <class Ty> Ty AbcdSquared_t<Ty>::operator()(Time t) const {
    return abcd_->covariance(t, T_, S_);
}
}

#endif
