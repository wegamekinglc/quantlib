/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
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

/*! \file swaptionvolcube2.hpp
    \brief Swaption volatility cube, fit-later-interpolate-early approach
*/

#ifndef quantlib_swaption_volcube_fit_later_interpolate_early_h
#define quantlib_swaption_volcube_fit_later_interpolate_early_h

#include <ql/termstructures/volatility/swaption/swaptionvolcube.hpp>
#include <ql/math/interpolations/interpolation2d.hpp>
#include <ql/termstructures/volatility/interpolatedsmilesection.hpp>
#include <ql/math/interpolations/bilinearinterpolation.hpp>
#include <ql/math/rounding.hpp>

namespace QuantLib {

template <class T>
class SwaptionVolCube2_t : public SwaptionVolatilityCube_t<T> {
  public:
    /*! The swaption vol cube is made up of ordered swaption vol surface
        layers, each layer referring to a swap index of a given length
        (in years), all indexes belonging to the same family. In order
        to identify the family (and its market conventions) an index of
        whatever length from that family must be passed in as
        swapIndexBase.

        Often for short swap length the swap index family is different,
        e.g. the EUR case: swap vs 6M Euribor is used for length>1Y,
        while swap vs 3M Euribor is used for the 1Y length. The
        shortSwapIndexBase is used to identify this second family.
  */
    SwaptionVolCube2_t(
        const Handle<SwaptionVolatilityStructure> &atmVolStructure,
        const std::vector<Period> &optionTenors,
        const std::vector<Period> &swapTenors,
        const std::vector<T> &strikeSpreads,
        const std::vector<std::vector<Handle<Quote_t<T> > > > &volSpreads,
        const boost::shared_ptr<SwapIndex> &swapIndexBase,
        const boost::shared_ptr<SwapIndex> &shortSwapIndexBase,
        bool vegaWeightedSmileFit);
    //! \name LazyObject interface
    //@{
    void performCalculations() const;
    //@}
    //! \name SwaptionVolatilityCube inspectors
    //@{
    const Matrix_t<T> &volSpreads(Size i) const { return volSpreadsMatrix_[i]; }
    boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(const Date &optionDate, const Period &swapTenor) const;
    boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(Time optionTime, Time swapLength) const;
    //@}
  private:
    mutable std::vector<Interpolation2D_t<T> > volSpreadsInterpolator_;
    mutable std::vector<Matrix_t<T> > volSpreadsMatrix_;
};

typedef SwaptionVolCube2_t<Real> SwaptionVolCube2;

// implementation

template <class T>
SwaptionVolCube2_t<T>::SwaptionVolCube2_t(
    const Handle<SwaptionVolatilityStructure> &atmVolStructure,
    const std::vector<Period> &optionTenors,
    const std::vector<Period> &swapTenors, const std::vector<T> &strikeSpreads,
    const std::vector<std::vector<Handle<Quote_t<T> > > > &volSpreads,
    const boost::shared_ptr<SwapIndex> &swapIndexBase,
    const boost::shared_ptr<SwapIndex> &shortSwapIndexBase,
    bool vegaWeightedSmileFit)
    : SwaptionVolatilityCube(atmVolStructure, optionTenors, swapTenors,
                             strikeSpreads, volSpreads, swapIndexBase,
                             shortSwapIndexBase, vegaWeightedSmileFit),
      volSpreadsInterpolator_(this->nStrikes_),
      volSpreadsMatrix_(this->nStrikes_, Matrix_t<T>(optionTenors.size(),
                                                     swapTenors.size(), 0.0)) {}

template <class T> void SwaptionVolCube2_t<T>::performCalculations() const {

    SwaptionVolatilityDiscrete::performCalculations();
    //! set volSpreadsMatrix_ by volSpreads_ quotes
    for (Size i = 0; i < this->nStrikes_; i++)
        for (Size j = 0; j < this->nOptionTenors_; j++)
            for (Size k = 0; k < this->nSwapTenors_; k++) {
                volSpreadsMatrix_[i][j][k] =
                    this->volSpreads_[j * this->nSwapTenors_ + k][i]->value();
            }
    //! create volSpreadsInterpolator_
    for (Size i = 0; i < this->nStrikes_; i++) {
        volSpreadsInterpolator_[i] = BilinearInterpolation(
            this->swapLengths_.begin(), this->swapLengths_.end(),
            this->optionTimes_.begin(), this->optionTimes_.end(),
            volSpreadsMatrix_[i]);
        volSpreadsInterpolator_[i].enableExtrapolation();
    }
}

template <class T>
boost::shared_ptr<SmileSection_t<T> >
SwaptionVolCube2_t<T>::smileSectionImpl(Time optionTime,
                                        Time swapLength) const {

    this->calculate();
    Date optionDate =
        Date(static_cast<BigInteger>(this->optionInterpolator_(optionTime)));
    Rounding rounder(0);
    Period swapTenor(static_cast<Integer>(rounder(swapLength * 12.0)), Months);
    return smileSectionImpl(optionDate, swapTenor);
}

template <class T>
boost::shared_ptr<SmileSection_t<T> >
SwaptionVolCube2_t<T>::smileSectionImpl(const Date &optionDate,
                                        const Period &swapTenor) const {
    this->calculate();
    T atmForward = this->atmStrike(optionDate, swapTenor);
    T atmVol = this->atmVol_->volatility(optionDate, swapTenor, atmForward);
    Time optionTime = this->timeFromReference(optionDate);
    T exerciseTimeSqrt = QLFCT::sqrt(optionTime);
    std::vector<T> strikes, stdDevs;
    strikes.reserve(this->nStrikes_);
    stdDevs.reserve(this->nStrikes_);
    Time length = this->swapLength(swapTenor);
    for (Size i = 0; i < this->nStrikes_; ++i) {
        strikes.push_back(atmForward + this->strikeSpreads_[i]);
        stdDevs.push_back(
            exerciseTimeSqrt *
            (atmVol + volSpreadsInterpolator_[i](length, optionTime)));
    }
    return boost::shared_ptr<SmileSection_t<T> >(
        new InterpolatedSmileSection_t<Linear, T>(optionTime, strikes, stdDevs,
                                                  atmForward));
}
}

#endif
