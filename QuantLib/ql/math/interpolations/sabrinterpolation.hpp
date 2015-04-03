/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
 Copyright (C) 2007 Marco Bianchetti
 Copyright (C) 2007 François du Vignaud
 Copyright (C) 2007 Giorgio Facchinetti
 Copyright (C) 2006 Mario Pucci
 Copyright (C) 2006 StatPro Italia srl
 Copyright (C) 2014, 2015 Peter Caspers

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

/*! \file sabrinterpolation.hpp
    \brief SABR interpolation interpolation between discrete points
*/

#ifndef quantlib_sabr_interpolation_hpp
#define quantlib_sabr_interpolation_hpp

#include <ql/math/interpolations/xabrinterpolation.hpp>
#include <ql/termstructures/volatility/sabr.hpp>

#include <boost/make_shared.hpp>
#include <boost/assign/list_of.hpp>

namespace QuantLib {

namespace detail {

template <class T> class SABRWrapper_t {
  public:
    SABRWrapper_t(const Time t, const T &forward, const std::vector<T> &params)
        : t_(t), forward_(forward), params_(params) {
        validateSabrParameters(params[0], params[1], params[2], params[3]);
    }
    T volatility(const T x) {
        return sabrVolatility(x, forward_, t_, params_[0], params_[1],
                              params_[2], params_[3]);
    }

  private:
    const T t_, &forward_;
    const std::vector<T> &params_;
};

template <class T> struct SABRSpecs_t {
    Size dimension() { return 4; }
    void defaultValues(std::vector<T> &params, std::vector<bool> &,
                       const T &forward, const T expiryTIme) {
        if (params[1] == Null<T>())
            params[1] = 0.5;
        if (params[0] == Null<T>())
            // adapt alpha to beta level
            params[0] =
                0.2 * (params[1] < 0.9999 ? QLFCT::pow(forward, 1.0 - params[1])
                                          : 1.0);
        if (params[2] == Null<T>())
            params[2] = QLFCT::sqrt(0.4);
        if (params[3] == Null<T>())
            params[3] = 0.0;
    }
    void guess(Array_t<T> &values, const std::vector<bool> &paramIsFixed,
               const T &forward, const T expiryTime, const std::vector<T> &r) {
        Size j = 0;
        if (!paramIsFixed[1])
            values[1] = (1.0 - 2E-6) * r[j++] + 1E-6;
        if (!paramIsFixed[0]) {
            values[0] = (1.0 - 2E-6) * r[j++] + 1E-6; // lognormal vol guess
            // adapt this to beta level
            if (values[1] < 0.999)
                values[0] *= QLFCT::pow(forward, 1.0 - values[1]);
        }
        if (!paramIsFixed[2])
            values[2] = 1.5 * r[j++] + 1E-6;
        if (!paramIsFixed[3])
            values[3] = (2.0 * r[j++] - 1.0) * (1.0 - 1E-6);
    }
    T eps1() { return .0000001; }
    T eps2() { return .9999; }
    T dilationFactor() { return 0.001; }
    Array_t<T> inverse(const Array_t<T> &y, const std::vector<bool> &,
                       const std::vector<T> &, const Real) {
        Array_t<T> x(4);
        x[0] = y[0] < 25.0 + eps1() ? QLFCT::sqrt(y[0] - eps1())
                                    : (y[0] - eps1() + 25.0) / 10.0;
        // y_[1] = std::tan(M_PI*(x[1] - 0.5))/dilationFactor();
        x[1] = QLFCT::sqrt(-QLFCT::log(y[1]));
        x[2] = y[2] < 25.0 + eps1() ? QLFCT::sqrt(y[2] - eps1())
                                    : (y[2] - eps1() + 25.0) / 10.0;
        x[3] = std::asin(y[3] / eps2());
        return x;
    }
    Array_t<T> direct(const Array_t<T> &x, const std::vector<bool> &,
                      const std::vector<T> &, const Real) {
        Array_t<T> y(4);
        y[0] = std::fabs(x[0]) < 5.0 ? x[0] * x[0] + eps1()
                                     : (10.0 * std::fabs(x[0]) - 25.0) + eps1();
        // y_[1] = std::atan(dilationFactor_*x[1])/M_PI + 0.5;
        y[1] = std::fabs(x[1]) < QLFCT::sqrt(-QLFCT::log(eps1()))
                   ? QLFCT::exp(-(x[1] * x[1]))
                   : eps1();
        y[2] = std::fabs(x[2]) < 5.0 ? x[2] * x[2] + eps1()
                                     : (10.0 * std::fabs(x[2]) - 25.0) + eps1();
        y[3] = std::fabs(x[3]) < 2.5 * M_PI
                   ? eps2() * std::sin(x[3])
                   : eps2() * (x[3] > 0.0 ? 1.0 : (-1.0));
        return y;
    }
    typedef SABRWrapper_t<T> type;
    boost::shared_ptr<type> instance(const Time t, const T &forward,
                                     const std::vector<T> &params) {
        return boost::make_shared<type>(t, forward, params);
    }
};
}

//! %SABR smile interpolation between discrete volatility points.
template <class T> class SABRInterpolation_t : public Interpolation_t<T> {
  public:
    template <class I1, class I2>
    SABRInterpolation_t(const I1 &xBegin, // x = strikes
                        const I1 &xEnd,
                        const I2 &yBegin, // y = volatilities
                        Time t,           // option expiry
                        const T &forward, T alpha, T beta, T nu, T rho,
                        bool alphaIsFixed, bool betaIsFixed, bool nuIsFixed,
                        bool rhoIsFixed, bool vegaWeighted = true,
                        const boost::shared_ptr<EndCriteria> &endCriteria =
                            boost::shared_ptr<EndCriteria>(),
                        const boost::shared_ptr<OptimizationMethod> &optMethod =
                            boost::shared_ptr<OptimizationMethod>(),
                        const T errorAccept = 0.0020,
                        const bool useMaxError = false,
                        const Size maxGuesses = 50) {

        this->impl_ = boost::shared_ptr<typename Interpolation_t<T>::Impl>(
            new detail::XABRInterpolationImpl_t<I1, I2, detail::SABRSpecs_t, T>(
                xBegin, xEnd, yBegin, t, forward,
                boost::assign::list_of(alpha)(beta)(nu)(rho),
                boost::assign::list_of(alphaIsFixed)(betaIsFixed)(nuIsFixed)(
                    rhoIsFixed),
                vegaWeighted, endCriteria, optMethod, errorAccept, useMaxError,
                maxGuesses));
        coeffs_ = boost::dynamic_pointer_cast<
            detail::XABRCoeffHolder_t<detail::SABRSpecs_t, T> >(this->impl_);
    }
    T expiry() const { return coeffs_->t_; }
    T forward() const { return coeffs_->forward_; }
    T alpha() const { return coeffs_->params_[0]; }
    T beta() const { return coeffs_->params_[1]; }
    T nu() const { return coeffs_->params_[2]; }
    T rho() const { return coeffs_->params_[3]; }
    T rmsError() const { return coeffs_->error_; }
    T maxError() const { return coeffs_->maxError_; }
    const std::vector<T> &interpolationWeights() const {
        return coeffs_->weights_;
    }
    EndCriteria::Type endCriteria() { return coeffs_->XABREndCriteria_; }

  private:
    boost::shared_ptr<detail::XABRCoeffHolder_t<detail::SABRSpecs_t, T> >
        coeffs_;
};

typedef SABRInterpolation_t<Real> SABRInterpolation;

//! %SABR interpolation factory and traits
template <class T> class SABR {
  public:
    SABR(Time t, T forward, T alpha, T beta, T nu, T rho, bool alphaIsFixed,
         bool betaIsFixed, bool nuIsFixed, bool rhoIsFixed,
         bool vegaWeighted = false,
         const boost::shared_ptr<EndCriteria> endCriteria =
             boost::shared_ptr<EndCriteria>(),
         const boost::shared_ptr<OptimizationMethod> optMethod =
             boost::shared_ptr<OptimizationMethod>(),
         const T errorAccept = 0.0020, const bool useMaxError = false,
         const Size maxGuesses = 50)
        : t_(t), forward_(forward), alpha_(alpha), beta_(beta), nu_(nu),
          rho_(rho), alphaIsFixed_(alphaIsFixed), betaIsFixed_(betaIsFixed),
          nuIsFixed_(nuIsFixed), rhoIsFixed_(rhoIsFixed),
          vegaWeighted_(vegaWeighted), endCriteria_(endCriteria),
          optMethod_(optMethod), errorAccept_(errorAccept),
          useMaxError_(useMaxError), maxGuesses_(maxGuesses) {}
    template <class I1, class I2>
    Interpolation_t<T> interpolate(const I1 &xBegin, const I1 &xEnd,
                                   const I2 &yBegin) const {
        return SABRInterpolation_t<T>(
            xBegin, xEnd, yBegin, t_, forward_, alpha_, beta_, nu_, rho_,
            alphaIsFixed_, betaIsFixed_, nuIsFixed_, rhoIsFixed_, vegaWeighted_,
            endCriteria_, optMethod_, errorAccept_, useMaxError_, maxGuesses_);
    }
    static const bool global = true;

  private:
    Time t_;
    T forward_;
    T alpha_, beta_, nu_, rho_;
    bool alphaIsFixed_, betaIsFixed_, nuIsFixed_, rhoIsFixed_;
    bool vegaWeighted_;
    const boost::shared_ptr<EndCriteria> endCriteria_;
    const boost::shared_ptr<OptimizationMethod> optMethod_;
    const T errorAccept_;
    const bool useMaxError_;
    const Size maxGuesses_;
};
}

#endif
