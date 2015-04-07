/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013, 2015 Peter Caspers

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

/*! \file gaussian1dswaptionengine.hpp
    \brief
*/

#ifndef quantlib_pricers_gaussian1d_swaption_hpp
#define quantlib_pricers_gaussian1d_swaption_hpp

#include <ql/instruments/swaption.hpp>
#include <ql/pricingengines/genericmodelengine.hpp>
#include <ql/experimental/models/gaussian1dmodel.hpp>

namespace QuantLib {

//! One factor model swaption engine
/*! \ingroup swaptionengines

    All fixed coupons with start date greater or equal to the respective
    option expiry are considered to be
    part of the exercise into right.

    \warning Cash settled swaptions are not supported
*/

template <class T>
class Gaussian1dSwaptionEngine_t
    : public GenericModelEngine<Gaussian1dModel_t<T>,
                                typename Swaption_t<T>::arguments,
                                typename Swaption_t<T>::results> {
  public:
    Gaussian1dSwaptionEngine_t(
        const boost::shared_ptr<Gaussian1dModel> &model,
        const int integrationPoints = 64, const T stddevs = 7.0,
        const bool extrapolatePayoff = true,
        const bool flatPayoffExtrapolation = false,
        const Handle<YieldTermStructure_t<T> > &discountCurve =
            Handle<YieldTermStructure_t<T> >())
        : GenericModelEngine<Gaussian1dModel, typename Swaption_t<T>::arguments,
                             typename Swaption_t<T>::results>(model),
          integrationPoints_(integrationPoints), stddevs_(stddevs),
          extrapolatePayoff_(extrapolatePayoff),
          flatPayoffExtrapolation_(flatPayoffExtrapolation),
          discountCurve_(discountCurve) {

        if (!discountCurve_.empty())
            this->registerWith(discountCurve_);
    }

    void calculate() const;

  private:
    const int integrationPoints_;
    const T stddevs_;
    const bool extrapolatePayoff_, flatPayoffExtrapolation_;
    const Handle<YieldTermStructure_t<T> > discountCurve_;
};

typedef Gaussian1dSwaptionEngine_t<Real> Gaussian1dSwaptionEngine;

// implementation

template <class T> void Gaussian1dSwaptionEngine_t<T>::calculate() const {

    QL_REQUIRE(this->arguments_.settlementType == Settlement::Physical,
               "cash-settled swaptions not yet implemented ...");

    Date settlement = this->model_->termStructure()->referenceDate();

    if (this->arguments_.exercise->dates().back() <=
        settlement) { // swaption is expired, possibly generated swap is not
                      // valued
        this->results_.value = 0.0;
        return;
    }

    int idx = static_cast<int>(this->arguments_.exercise->dates().size()) - 1;
    int minIdxAlive = static_cast<int>(
        std::upper_bound(this->arguments_.exercise->dates().begin(),
                         this->arguments_.exercise->dates().end(), settlement) -
        this->arguments_.exercise->dates().begin());

    VanillaSwap_t<T> swap = *this->arguments_.swap;
    Option::Type type = this->arguments_.type == VanillaSwap_t<T>::Payer
                            ? Option::Call
                            : Option::Put;
    Schedule fixedSchedule = swap.fixedSchedule();
    Schedule floatSchedule = swap.floatingSchedule();

    Array_t<T> npv0(2 * integrationPoints_ + 1, 0.0),
        npv1(2 * integrationPoints_ + 1, 0.0);
    Array_t<T> z = this->model_->yGrid(stddevs_, integrationPoints_);
    Array_t<T> p(z.size(), 0.0);

    Date expiry1 = Null<Date>(), expiry0;
    Time expiry1Time = Null<T>(), expiry0Time;

    do {

        if (idx == minIdxAlive - 1)
            expiry0 = settlement;
        else
            expiry0 = this->arguments_.exercise->dates()[idx];

        expiry0Time = QLFCT::max(
            this->model_->termStructure()->timeFromReference(expiry0), 0.0);

        Size j1 = std::upper_bound(fixedSchedule.dates().begin(),
                                   fixedSchedule.dates().end(), expiry0 - 1) -
                  fixedSchedule.dates().begin();
        Size k1 = std::upper_bound(floatSchedule.dates().begin(),
                                   floatSchedule.dates().end(), expiry0 - 1) -
                  floatSchedule.dates().begin();

// a lazy object is not thread safe, neither is the caching
// in gsrprocess. therefore we trigger computations here such
// that neither lazy object recalculation nor write access
// during caching occurs in the parallized loop below.
// this is known to work for the gsr and markov functional
// model implementations of Gaussian1dModel
#ifdef _OPENMP
        if (expiry0 > settlement) {
            for (Size l = k1; l < this->arguments_.floatingCoupons.size();
                 l++) {
                this->model_->forwardRate(this->arguments_.floatingFixingDates[l],
                                    expiry0, 0.0,
                                    this->arguments_.swap->iborIndex());
                this->model_->zerobond(this->arguments_.floatingPayDates[l], expiry0,
                                 0.0, discountCurve_);
            }
            for (Size l = j1; l < this->arguments_.fixedCoupons.size(); l++) {
                this->model_->zerobond(this->arguments_.fixedPayDates[l], expiry0,
                                 0.0, discountCurve_);
            }
            this->model_->numeraire(expiry0Time, 0.0, discountCurve_);
        }
#endif

#pragma omp parallel for default(shared)                                       \
    firstprivate(p) if (expiry0 > settlement)
        for (Size k = 0; k < (expiry0 > settlement ? npv0.size() : 1); k++) {

            T price = 0.0;
            if (expiry1Time != Null<T>()) {
                Array_t<T> yg = this->model_->yGrid(
                    stddevs_, integrationPoints_, expiry1Time, expiry0Time,
                    expiry0 > settlement ? z[k] : 0.0);
                CubicInterpolation payoff0(z.begin(), z.end(), npv1.begin(),
                                           CubicInterpolation::Spline, true,
                                           CubicInterpolation::Lagrange, 0.0,
                                           CubicInterpolation::Lagrange, 0.0);
                for (Size i = 0; i < yg.size(); i++) {
                    p[i] = payoff0(yg[i], true);
                }
                CubicInterpolation payoff1(z.begin(), z.end(), p.begin(),
                                           CubicInterpolation::Spline, true,
                                           CubicInterpolation::Lagrange, 0.0,
                                           CubicInterpolation::Lagrange, 0.0);
                for (Size i = 0; i < z.size() - 1; i++) {
                    price += this->model_->gaussianShiftedPolynomialIntegral(
                        0.0, payoff1.cCoefficients()[i],
                        payoff1.bCoefficients()[i], payoff1.aCoefficients()[i],
                        p[i], z[i], z[i], z[i + 1]);
                }
                if (extrapolatePayoff_) {
                    if (flatPayoffExtrapolation_) {
                        price += this->model_->gaussianShiftedPolynomialIntegral(
                            0.0, 0.0, 0.0, 0.0, p[z.size() - 2],
                            z[z.size() - 2], z[z.size() - 1], 100.0);
                        price += this->model_->gaussianShiftedPolynomialIntegral(
                            0.0, 0.0, 0.0, 0.0, p[0], z[0], -100.0, z[0]);
                    } else {
                        if (type == Option::Call)
                            price += this->model_->gaussianShiftedPolynomialIntegral(
                                0.0, payoff1.cCoefficients()[z.size() - 2],
                                payoff1.bCoefficients()[z.size() - 2],
                                payoff1.aCoefficients()[z.size() - 2],
                                p[z.size() - 2], z[z.size() - 2],
                                z[z.size() - 1], 100.0);
                        if (type == Option::Put)
                            price += this->model_->gaussianShiftedPolynomialIntegral(
                                0.0, payoff1.cCoefficients()[0],
                                payoff1.bCoefficients()[0],
                                payoff1.aCoefficients()[0], p[0], z[0], -100.0,
                                z[0]);
                    }
                }
            }

            npv0[k] = price;

            if (expiry0 > settlement) {
                T floatingLegNpv = 0.0;
                for (Size l = k1; l < this->arguments_.floatingCoupons.size();
                     l++) {
                    floatingLegNpv +=
                        this->arguments_.nominal *
                        this->arguments_.floatingAccrualTimes[l] *
                        (this->arguments_.floatingSpreads[l] +
                         this->model_->forwardRate(
                             this->arguments_.floatingFixingDates[l], expiry0,
                             z[k], this->arguments_.swap->iborIndex())) *
                        this->model_->zerobond(this->arguments_.floatingPayDates[l],
                                         expiry0, z[k], discountCurve_);
                }
                T fixedLegNpv = 0.0;
                for (Size l = j1; l < this->arguments_.fixedCoupons.size();
                     l++) {
                    fixedLegNpv +=
                        this->arguments_.fixedCoupons[l] *
                        this->model_->zerobond(this->arguments_.fixedPayDates[l],
                                         expiry0, z[k], discountCurve_);
                }
                npv0[k] =
                    QLFCT::max(npv0[k], (type == Option::Call ? 1.0 : -1.0) *
                                            (floatingLegNpv - fixedLegNpv) /
                                            this->model_->numeraire(expiry0Time, z[k],
                                                              discountCurve_));
            }
        }

        npv1.swap(npv0);
        expiry1 = expiry0;
        expiry1Time = expiry0Time;

    } while (--idx >= minIdxAlive - 1);

    this->results_.value =
        npv1[0] * this->model_->numeraire(0.0, 0.0, discountCurve_);
}

} // namespace QuantLib

#endif
