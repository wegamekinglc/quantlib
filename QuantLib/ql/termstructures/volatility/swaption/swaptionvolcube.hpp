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

/*! \file swaptionvolcube.hpp
    \brief Swaption volatility cube
*/

#ifndef quantlib_swaption_volatility_cube_h
#define quantlib_swaption_volatility_cube_h

#include <ql/termstructures/volatility/swaption/swaptionvoldiscrete.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>
#include <ql/quote.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/indexes/swapindex.hpp>

namespace QuantLib {

//! swaption-volatility cube
/*! \warning this class is not finalized and its interface might
             change in subsequent releases.
*/
template <class T>
class SwaptionVolatilityCube_t : public SwaptionVolatilityDiscrete_t<T> {
  public:
    SwaptionVolatilityCube_t(
        const Handle<SwaptionVolatilityStructure_t<T> > &atmVolStructure,
        const std::vector<Period> &optionTenors,
        const std::vector<Period> &swapTenors,
        const std::vector<T> &strikeSpreads,
        const std::vector<std::vector<Handle<Quote_t<T> > > > &volSpreads,
        const boost::shared_ptr<SwapIndex_t<T> > &swapIndexBase,
        const boost::shared_ptr<SwapIndex_t<T> > &shortSwapIndexBase,
        bool vegaWeightedSmileFit);
    //! \name TermStructure interface
    //@{
    DayCounter dayCounter() const { return atmVol_->dayCounter(); }
    Date maxDate() const { return atmVol_->maxDate(); }
    Time maxTime() const { return atmVol_->maxTime(); }
    const Date &referenceDate() const { return atmVol_->referenceDate(); }
    Calendar calendar() const { return atmVol_->calendar(); }
    Natural settlementDays() const { return atmVol_->settlementDays(); }
    //! \name VolatilityTermStructure interface
    //@{
    T minStrike() const { return 0.0; }
    T maxStrike() const { return 1.0; }
    //@}
    //! \name SwaptionVolatilityStructure interface
    //@{
    const Period &maxSwapTenor() const { return atmVol_->maxSwapTenor(); }
    //@}
    //! \name Other inspectors
    //@{
    T atmStrike(const Date &optionDate, const Period &swapTenor) const;
    T atmStrike(const Period &optionTenor, const Period &swapTenor) const {
        Date optionDate = this->optionDateFromTenor(optionTenor);
        return atmStrike(optionDate, swapTenor);
    }
    //@}
  protected:
    void registerWithVolatilitySpread();
    T volatilityImpl(Time optionTime, Time swapLength, T strike) const;
    T volatilityImpl(const Date &optionDate, const Period &swapTenor,
                     T strike) const;
    Handle<SwaptionVolatilityStructure_t<T> > atmVol_;
    Size nStrikes_;
    std::vector<T> strikeSpreads_;
    mutable std::vector<T> localStrikes_;
    mutable std::vector<T> localSmile_;
    std::vector<std::vector<Handle<Quote_t<T> > > > volSpreads_;
    boost::shared_ptr<SwapIndex_t<T> > swapIndexBase_, shortSwapIndexBase_;
    bool vegaWeightedSmileFit_;
};

typedef SwaptionVolatilityCube_t<Real> SwaptionVolatilityCube;

// inline

template <class T>
inline T SwaptionVolatilityCube_t<T>::volatilityImpl(Time optionTime,
                                                     Time swapLength,
                                                     T strike) const {
    return this->smileSectionImpl(optionTime, swapLength)->volatility(strike);
}

template <class T>
inline T SwaptionVolatilityCube_t<T>::volatilityImpl(const Date &optionDate,
                                                     const Period &swapTenor,
                                                     T strike) const {
    return this->smileSectionImpl(optionDate, swapTenor)->volatility(strike);
}

// implementation

template <class T>
SwaptionVolatilityCube_t<T>::SwaptionVolatilityCube_t(
    const Handle<SwaptionVolatilityStructure_t<T> > &atmVol,
    const std::vector<Period> &optionTenors,
    const std::vector<Period> &swapTenors, const std::vector<T> &strikeSpreads,
    const std::vector<std::vector<Handle<Quote_t<T> > > > &volSpreads,
    const boost::shared_ptr<SwapIndex_t<T> > &swapIndexBase,
    const boost::shared_ptr<SwapIndex_t<T> > &shortSwapIndexBase,
    bool vegaWeightedSmileFit)
    : SwaptionVolatilityDiscrete_t<T>(
          optionTenors, swapTenors, 0, atmVol->calendar(),
          atmVol->businessDayConvention(), atmVol->dayCounter()),
      atmVol_(atmVol), nStrikes_(strikeSpreads.size()),
      strikeSpreads_(strikeSpreads), localStrikes_(nStrikes_),
      localSmile_(nStrikes_), volSpreads_(volSpreads),
      swapIndexBase_(swapIndexBase), shortSwapIndexBase_(shortSwapIndexBase),
      vegaWeightedSmileFit_(vegaWeightedSmileFit) {
    QL_REQUIRE(!atmVol.empty(), "atm vol handle not linked to anything");

    QL_REQUIRE(nStrikes_ >= 1, "too few strikes ("
                                   << nStrikes_
                                   << ")"); // quick fix (merge master !)
    for (Size i = 1; i < nStrikes_; ++i)
        QL_REQUIRE(strikeSpreads_[i - 1] < strikeSpreads_[i],
                   "non increasing strike spreads: "
                       << io::ordinal(i) << " is " << strikeSpreads_[i - 1]
                       << ", " << io::ordinal(i + 1) << " is "
                       << strikeSpreads_[i]);

    QL_REQUIRE(!volSpreads_.empty(), "empty vol spreads matrix");

    QL_REQUIRE(this->nOptionTenors_ * this->nSwapTenors_ == volSpreads_.size(),
               "mismatch between number of option tenors * swap tenors ("
                   << this->nOptionTenors_ * this->nSwapTenors_
                   << ") and number of rows (" << volSpreads_.size() << ")");

    for (Size i = 0; i < volSpreads_.size(); i++)
        QL_REQUIRE(this->nStrikes_ == this->volSpreads_[i].size(),
                   "mismatch between number of strikes ("
                       << this->nStrikes_ << ") and number of columns ("
                       << volSpreads_[i].size() << ") in the "
                       << io::ordinal(i + 1) << " row");

    this->registerWith(atmVol_);
    atmVol_->enableExtrapolation();

    this->registerWith(swapIndexBase_);
    this->registerWith(shortSwapIndexBase_);

    QL_REQUIRE(shortSwapIndexBase_->tenor() < swapIndexBase_->tenor(),
               "short index tenor (" << shortSwapIndexBase_->tenor()
                                     << ") is not less than index tenor ("
                                     << swapIndexBase_->tenor() << ")");

    this->registerWithVolatilitySpread();
    this->registerWith(Settings::instance().evaluationDate());
    this->evaluationDate_ = Settings::instance().evaluationDate();
}

template <class T>
void SwaptionVolatilityCube_t<T>::registerWithVolatilitySpread() {
    for (Size i = 0; i < nStrikes_; i++)
        for (Size j = 0; j < this->nOptionTenors_; j++)
            for (Size k = 0; k < this->nSwapTenors_; k++)
                this->registerWith(volSpreads_[j * this->nSwapTenors_ + k][i]);
}

template <class T>
T SwaptionVolatilityCube_t<T>::atmStrike(const Date &optionD,
                                         const Period &swapTenor) const {

    // FIXME use a familyName-based index factory
    if (swapTenor > shortSwapIndexBase_->tenor()) {
        if (swapIndexBase_->exogenousDiscount()) {
            return SwapIndex_t<T>(swapIndexBase_->familyName(), swapTenor,
                                  swapIndexBase_->fixingDays(),
                                  swapIndexBase_->currency(),
                                  swapIndexBase_->fixingCalendar(),
                                  swapIndexBase_->fixedLegTenor(),
                                  swapIndexBase_->fixedLegConvention(),
                                  swapIndexBase_->dayCounter(),
                                  swapIndexBase_->iborIndex(),
                                  swapIndexBase_->discountingTermStructure())
                .fixing(optionD);
        } else {
            return SwapIndex_t<T>(swapIndexBase_->familyName(), swapTenor,
                                  swapIndexBase_->fixingDays(),
                                  swapIndexBase_->currency(),
                                  swapIndexBase_->fixingCalendar(),
                                  swapIndexBase_->fixedLegTenor(),
                                  swapIndexBase_->fixedLegConvention(),
                                  swapIndexBase_->dayCounter(),
                                  swapIndexBase_->iborIndex())
                .fixing(optionD);
        }
    } else {
        if (shortSwapIndexBase_->exogenousDiscount()) {
            return SwapIndex_t<T>(
                       shortSwapIndexBase_->familyName(), swapTenor,
                       shortSwapIndexBase_->fixingDays(),
                       shortSwapIndexBase_->currency(),
                       shortSwapIndexBase_->fixingCalendar(),
                       shortSwapIndexBase_->fixedLegTenor(),
                       shortSwapIndexBase_->fixedLegConvention(),
                       shortSwapIndexBase_->dayCounter(),
                       shortSwapIndexBase_->iborIndex(),
                       shortSwapIndexBase_->discountingTermStructure())
                .fixing(optionD);
        } else {
            return SwapIndex_t<T>(shortSwapIndexBase_->familyName(), swapTenor,
                                  shortSwapIndexBase_->fixingDays(),
                                  shortSwapIndexBase_->currency(),
                                  shortSwapIndexBase_->fixingCalendar(),
                                  shortSwapIndexBase_->fixedLegTenor(),
                                  shortSwapIndexBase_->fixedLegConvention(),
                                  shortSwapIndexBase_->dayCounter(),
                                  shortSwapIndexBase_->iborIndex())
                .fixing(optionD);
        }
    }
}
}

#endif
