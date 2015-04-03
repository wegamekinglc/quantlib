/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2007 Cristina Duminuco
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

/*! \file sabrinterpolatedsmilesection.hpp
    \brief Interpolated smile section class
*/

#ifndef quantlib_sabr_interpolated_smile_section_hpp
#define quantlib_sabr_interpolated_smile_section_hpp

#include <ql/handle.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>
#include <ql/math/interpolations/sabrinterpolation.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/settings.hpp>
#include <ql/quotes/simplequote.hpp>

namespace QuantLib {

template <class T>
class SabrInterpolatedSmileSection_t : public SmileSection_t<T>,
                                       public LazyObject {
  public:
    //! \name Constructors
    //@{
    //! all market data are quotes
    SabrInterpolatedSmileSection_t(
        const Date &optionDate, const Handle<Quote_t<T> > &forward,
        const std::vector<T> &strikes, bool hasFloatingStrikes,
        const Handle<Quote_t<T> > &atmVolatility,
        const std::vector<Handle<Quote_t<T> > > &volHandles, T alpha, T beta,
        T nu, T rho, bool isAlphaFixed = false, bool isBetaFixed = false,
        bool isNuFixed = false, bool isRhoFixed = false,
        bool vegaWeighted = true,
        const boost::shared_ptr<EndCriteria> &endCriteria =
            boost::shared_ptr<EndCriteria>(),
        const boost::shared_ptr<OptimizationMethod> &method =
            boost::shared_ptr<OptimizationMethod>(),
        const DayCounter &dc = Actual365Fixed());
    //! no quotes
    SabrInterpolatedSmileSection_t(
        const Date &optionDate, const T &forward, const std::vector<T> &strikes,
        bool hasFloatingStrikes, const T &atmVolatility,
        const std::vector<T> &vols, T alpha, T beta, T nu, T rho,
        bool isAlphaFixed = false, bool isBetaFixed = false,
        bool isNuFixed = false, bool isRhoFixed = false,
        bool vegaWeighted = true,
        const boost::shared_ptr<EndCriteria> &endCriteria =
            boost::shared_ptr<EndCriteria>(),
        const boost::shared_ptr<OptimizationMethod> &method =
            boost::shared_ptr<OptimizationMethod>(),
        const DayCounter &dc = Actual365Fixed());
    //@}
    //! \name LazyObject interface
    //@{
    virtual void performCalculations() const;
    virtual void update();
    //@}
    //! \name SmileSection interface
    //@{
    T minStrike() const;
    T maxStrike() const;
    T atmLevel() const;
    //@}
    T varianceImpl(T strike) const;
    T volatilityImpl(T strike) const;
    //! \name Inspectors
    //@{
    T alpha() const;
    T beta() const;
    T nu() const;
    T rho() const;
    T rmsError() const;
    T maxError() const;
    EndCriteria::Type endCriteria() const;
    //@}

  protected:
    //! Creates the mutable SABRInterpolation
    void createInterpolation() const;
    mutable boost::shared_ptr<SABRInterpolation_t<T>> sabrInterpolation_;

    //! Market data
    const Handle<Quote_t<T> > forward_;
    const Handle<Quote_t<T> > atmVolatility_;
    std::vector<Handle<Quote_t<T> > > volHandles_;
    mutable std::vector<T> strikes_;
    //! Only strikes corresponding to valid market data
    mutable std::vector<T> actualStrikes_;
    bool hasFloatingStrikes_;

    mutable T forwardValue_;
    mutable std::vector<T> vols_;
    //! Sabr parameters
    T alpha_, beta_, nu_, rho_;
    //! Sabr interpolation settings
    bool isAlphaFixed_, isBetaFixed_, isNuFixed_, isRhoFixed_;
    bool vegaWeighted_;
    const boost::shared_ptr<EndCriteria> endCriteria_;
    const boost::shared_ptr<OptimizationMethod> method_;

    mutable Date evaluationDate_;
};

typedef SabrInterpolatedSmileSection_t<Real> SabrInterpolatedSmileSection;

template <class T> inline void SabrInterpolatedSmileSection_t<T>::update() {
    LazyObject::update();
    SmileSection::update();
}

template <class T>
inline T SabrInterpolatedSmileSection_t<T>::volatilityImpl(T strike) const {
    calculate();
    return (*sabrInterpolation_)(strike, true);
}

template <class T> inline T SabrInterpolatedSmileSection_t<T>::alpha() const {
    calculate();
    return sabrInterpolation_->alpha();
}

template <class T> inline T SabrInterpolatedSmileSection_t<T>::beta() const {
    calculate();
    return sabrInterpolation_->beta();
}

template <class T> inline T SabrInterpolatedSmileSection_t<T>::nu() const {
    calculate();
    return sabrInterpolation_->nu();
}

template <class T> inline T SabrInterpolatedSmileSection_t<T>::rho() const {
    calculate();
    return sabrInterpolation_->rho();
}

template <class T>
inline T SabrInterpolatedSmileSection_t<T>::rmsError() const {
    calculate();
    return sabrInterpolation_->rmsError();
}

template <class T>
inline T SabrInterpolatedSmileSection_t<T>::maxError() const {
    calculate();
    return sabrInterpolation_->maxError();
}

template <class T>
inline EndCriteria::Type
SabrInterpolatedSmileSection_t<T>::endCriteria() const {
    calculate();
    return sabrInterpolation_->endCriteria();
}

template <class T>
inline T SabrInterpolatedSmileSection_t<T>::minStrike() const {
    calculate();
    return actualStrikes_.front();
}

template <class T>
inline T SabrInterpolatedSmileSection_t<T>::maxStrike() const {
    calculate();
    return actualStrikes_.back();
}

template <class T>
inline T SabrInterpolatedSmileSection_t<T>::atmLevel() const {
    calculate();
    return forwardValue_;
}

// implementation

template <class T>
SabrInterpolatedSmileSection_t<T>::SabrInterpolatedSmileSection_t(
    const Date &optionDate, const Handle<Quote_t<T> > &forward,
    const std::vector<T> &strikes, bool hasFloatingStrikes,
    const Handle<Quote_t<T> > &atmVolatility,
    const std::vector<Handle<Quote_t<T> > > &volHandles, T alpha, T beta, T nu,
    T rho, bool isAlphaFixed, bool isBetaFixed, bool isNuFixed, bool isRhoFixed,
    bool vegaWeighted, const boost::shared_ptr<EndCriteria> &endCriteria,
    const boost::shared_ptr<OptimizationMethod> &method, const DayCounter &dc)
    : SmileSection_t<T>(optionDate, dc), forward_(forward),
      atmVolatility_(atmVolatility), volHandles_(volHandles), strikes_(strikes),
      actualStrikes_(strikes), hasFloatingStrikes_(hasFloatingStrikes),
      vols_(volHandles.size()), alpha_(alpha), beta_(beta), nu_(nu), rho_(rho),
      isAlphaFixed_(isAlphaFixed), isBetaFixed_(isBetaFixed),
      isNuFixed_(isNuFixed), isRhoFixed_(isRhoFixed),
      vegaWeighted_(vegaWeighted), endCriteria_(endCriteria), method_(method),
      evaluationDate_(Settings::instance().evaluationDate()) {

    LazyObject::registerWith(forward_);
    LazyObject::registerWith(atmVolatility_);
    for (Size i = 0; i < volHandles_.size(); ++i)
        LazyObject::registerWith(volHandles_[i]);
}

template <class T>
SabrInterpolatedSmileSection_t<T>::SabrInterpolatedSmileSection_t(
    const Date &optionDate, const T &forward, const std::vector<T> &strikes,
    bool hasFloatingStrikes, const T &atmVolatility,
    const std::vector<T> &volHandles, T alpha, T beta, T nu, T rho,
    bool isAlphaFixed, bool isBetaFixed, bool isNuFixed, bool isRhoFixed,
    bool vegaWeighted, const boost::shared_ptr<EndCriteria> &endCriteria,
    const boost::shared_ptr<OptimizationMethod> &method, const DayCounter &dc)
    : SmileSection_t<T>(optionDate, dc),
      forward_(Handle<Quote_t<T> >(
          boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(forward)))),
      atmVolatility_(Handle<Quote_t<T> >(
          boost::shared_ptr<Quote_t<T> >(new SimpleQuote_t<T>(atmVolatility)))),
      volHandles_(volHandles.size()), strikes_(strikes),
      actualStrikes_(strikes), hasFloatingStrikes_(hasFloatingStrikes),
      vols_(volHandles.size()), alpha_(alpha), beta_(beta), nu_(nu), rho_(rho),
      isAlphaFixed_(isAlphaFixed), isBetaFixed_(isBetaFixed),
      isNuFixed_(isNuFixed), isRhoFixed_(isRhoFixed),
      vegaWeighted_(vegaWeighted), endCriteria_(endCriteria), method_(method),
      evaluationDate_(Settings::instance().evaluationDate()) {

    for (Size i = 0; i < volHandles_.size(); ++i)
        volHandles_[i] = Handle<Quote_t<T> >(boost::shared_ptr<Quote_t<T> >(
            new SimpleQuote_t<T>(volHandles[i])));
}

template <class T>
void SabrInterpolatedSmileSection_t<T>::createInterpolation() const {
    boost::shared_ptr<SABRInterpolation_t<T> > tmp(new SABRInterpolation_t<T>(
        actualStrikes_.begin(), actualStrikes_.end(), vols_.begin(),
        this->exerciseTime(), forwardValue_, alpha_, beta_, nu_, rho_, isAlphaFixed_,
        isBetaFixed_, isNuFixed_, isRhoFixed_, vegaWeighted_, endCriteria_,
        method_));
    swap(tmp, sabrInterpolation_);
}

template <class T>
void SabrInterpolatedSmileSection_t<T>::performCalculations() const {
    forwardValue_ = forward_->value();
    vols_.clear();
    actualStrikes_.clear();
    // we populate the volatilities, skipping the invalid ones
    for (Size i = 0; i < volHandles_.size(); ++i) {
        if (volHandles_[i]->isValid()) {
            if (hasFloatingStrikes_) {
                actualStrikes_.push_back(forwardValue_ + strikes_[i]);
                vols_.push_back(atmVolatility_->value() +
                                volHandles_[i]->value());
            } else {
                actualStrikes_.push_back(strikes_[i]);
                vols_.push_back(volHandles_[i]->value());
            }
        }
    }
    // we are recreating the sabrinterpolation object unconditionnaly to
    // avoid iterator invalidation
    createInterpolation();
    sabrInterpolation_->update();
}

template <class T>
T SabrInterpolatedSmileSection_t<T>::varianceImpl(T strike) const {
    calculate();
    T v = (*sabrInterpolation_)(strike, true);
    return v * v * this->exerciseTime();
}
}

#endif
