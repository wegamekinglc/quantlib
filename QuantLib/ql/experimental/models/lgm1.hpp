/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

/*! \file lgm1.hpp
    \brief LGM model with piecewise alpha and constant kappa
*/

#ifndef quantlib_lgm1_hpp
#define quantlib_lgm1_hpp

#include <ql/experimental/models/lgm.hpp>
#include <ql/experimental/models/lgmpiecewisealphaconstantkappa.hpp>
#include <ql/models/model.hpp>

namespace QuantLib {

class Lgm1 : public Lgm<detail::LgmPiecewiseAlphaConstantKappa>,
             public CalibratedModel {
  public:
    // fixed model data
    Lgm1(const Handle<YieldTermStructure> &yts,
         const std::vector<Date> &volstepdates, const std::vector<Real> &alphas,
         const Real &kappa);

    // floating model data
    Lgm1(const Handle<YieldTermStructure> &yts,
         const std::vector<Date> &volstepdates,
         const std::vector<Handle<Quote> > &alphas, const Handle<Quote> &kappa);

    const Array &alpha() const { return alpha_.params(); };
    const Real kappa() const { return kappa_.params()[0]; };

  protected:
    void generateArguments() { notifyObservers(); }
    void update() { LazyObject::update(); }
    void performCalculations() const {
        Lgm::performCalculations();
        updateTimes();
        parametrization_.update();
    }

  private:
    void updateTimes() const;
    void updateAlpha();
    void updateKappa();
    void initialize();

    const std::vector<Date> volstepdates_;
    mutable std::vector<Real>
        volsteptimes_;                // used for parameters, parametrization_
    mutable Array volsteptimesArray_; // used for state process

    std::vector<Handle<Quote> > alphaQuotes_; // floating model data
    Handle<Quote> kappaQuote_;

    Parameter &alpha_, &kappa_; // calibrated model parameters

    const detail::LgmPiecewiseAlphaConstantKappa parametrization_;

    struct AlphaObserver : public Observer {
        AlphaObserver(Lgm1 *p) : p_(p) {}
        void update() { p_->updateAlpha(); }
        Lgm1 *p_;
    };

    struct KappaObserver : public Observer {
        KappaObserver(Lgm1 *p) : p_(p) {}
        void update() { p_->updateKappa(); }
        Lgm1 *p_;
    };

    boost::shared_ptr<AlphaObserver> alphaObserver_;
    boost::shared_ptr<KappaObserver> kappaObserver_;
};

} // namespace QuantLib

#endif
