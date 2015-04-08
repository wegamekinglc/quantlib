/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Giorgio Facchinetti
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

/*! \file swaptionvolcube1.hpp
    \brief Swaption volatility cube, fit-early-interpolate-later approach
           The provided types are
           SwaptionVolCube1 using the classic Hagan 2002 Sabr formula
           SwaptionVolCube1a using the No Arbitrage Sabr model (Doust)
*/

#ifndef quantlib_swaption_volcube_fit_early_interpolate_later_h
#define quantlib_swaption_volcube_fit_early_interpolate_later_h

#include <ql/termstructures/volatility/swaption/swaptionvolcube.hpp>
#include <ql/termstructures/volatility/sabrsmilesection.hpp>
#include <ql/math/matrix.hpp>
#include <ql/math/interpolations/sabrinterpolation.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/math/interpolations/flatextrapolation2d.hpp>
#include <ql/math/interpolations/backwardflatlinearinterpolation.hpp>
#include <ql/math/interpolations/bilinearinterpolation.hpp>
#include <ql/experimental/volatility/noarbsabrinterpolation.hpp>
#include <ql/quote.hpp>

#include <boost/make_shared.hpp>

#ifndef SWAPTIONVOLCUBE_VEGAWEIGHTED_TOL
#define SWAPTIONVOLCUBE_VEGAWEIGHTED_TOL 15.0e-4
#endif
#ifndef SWAPTIONVOLCUBE_TOL
#define SWAPTIONVOLCUBE_TOL 100.0e-4
#endif

namespace QuantLib {

template <class> class Interpolation2D_t;
template <class> class EndCriteria_t;
template <class> class OptimizationMethod_t;

template <template <class> class Model, class T>
class SwaptionVolCube1x_t : public SwaptionVolatilityCube {
    class Cube {
      public:
        Cube() {}
        Cube(const std::vector<Date> &optionDates,
             const std::vector<Period> &swapTenors,
             const std::vector<Time> &optionTimes,
             const std::vector<Time> &swapLengths, Size nLayers,
             bool extrapolation = true, bool backwardFlat = false);
        Cube &operator=(const Cube &o);
        Cube(const Cube &);
        virtual ~Cube() {}
        void setElement(Size IndexOfLayer, Size IndexOfRow, Size IndexOfColumn,
                        T x);
        void setPoints(const std::vector<Matrix_t<T> > &x);
        void setPoint(const Date &optionDate, const Period &swapTenor,
                      const Time optionTime, const Time swapLengths,
                      const std::vector<T> &point);
        void setLayer(Size i, const Matrix_t<T> &x);
        void expandLayers(Size i, bool expandOptionTimes, Size j,
                          bool expandSwapLengths);
        const std::vector<Date> &optionDates() const { return optionDates_; }
        const std::vector<Period> &swapTenors() const { return swapTenors_; }
        const std::vector<Time> &optionTimes() const;
        const std::vector<Time> &swapLengths() const;
        const std::vector<Matrix_t<T> > &points() const;
        std::vector<T> operator()(const Time optionTime,
                                  const Time swapLengths) const;
        void updateInterpolators() const;
        Matrix_t<T> browse() const;

      private:
        std::vector<Time> optionTimes_, swapLengths_;
        std::vector<Date> optionDates_;
        std::vector<Period> swapTenors_;
        Size nLayers_;
        std::vector<Matrix_t<T> > points_;
        mutable std::vector<Matrix_t<T> > transposedPoints_;
        bool extrapolation_;
        bool backwardFlat_;
        mutable std::vector<boost::shared_ptr<Interpolation2D_t<T> > >
            interpolators_;
    };

  public:
    SwaptionVolCube1x_t(
        const Handle<SwaptionVolatilityStructure_t<T>> &atmVolStructure,
        const std::vector<Period> &optionTenors,
        const std::vector<Period> &swapTenors,
        const std::vector<T> &strikeSpreads,
        const std::vector<std::vector<Handle<Quote_t<T> > > > &volSpreads,
        const boost::shared_ptr<SwapIndex_t<T> > &swapIndexBase,
        const boost::shared_ptr<SwapIndex_t<T> > &shortSwapIndexBase,
        bool vegaWeightedSmileFit,
        const std::vector<std::vector<Handle<Quote_t<T> > > > &parametersGuess,
        const std::vector<bool> &isParameterFixed, bool isAtmCalibrated,
        const boost::shared_ptr<EndCriteria_t<T>> &endCriteria =
            boost::shared_ptr<EndCriteria_t<T>>(),
        T maxErrorTolerance = Null<T>(),
        const boost::shared_ptr<OptimizationMethod_t<T>> &optMethod =
            boost::shared_ptr<OptimizationMethod_t<T>>(),
        const T errorAccept = Null<T>(), const bool useMaxError = false,
        const Size maxGuesses = 50, const bool backwardFlat = false,
        const T cutoffStrike = 0.0001);
    //! \name LazyObject interface
    //@{
    void performCalculations() const;
    //@}
    //! \name SwaptionVolatilityCube interface
    //@{
    boost::shared_ptr<SmileSection_t<T> >
    smileSectionImpl(Time optionTime, Time swapLength) const;
    //@}
    //! \name Other inspectors
    //@{
    const Matrix_t<T> &marketVolCube(Size i) const {
        return marketVolCube_.points()[i];
    }
    Matrix_t<T> sparseSabrParameters() const;
    Matrix_t<T> denseSabrParameters() const;
    Matrix_t<T> marketVolCube() const;
    Matrix_t<T> volCubeAtmCalibrated() const;
    //@}
    void sabrCalibrationSection(const Cube &marketVolCube, Cube &parametersCube,
                                const Period &swapTenor) const;
    void recalibration(T beta, const Period &swapTenor);
    void recalibration(const std::vector<T> &beta, const Period &swapTenor);
    void recalibration(const std::vector<Period> &swapLengths,
                       const std::vector<T> &beta, const Period &swapTenor);
    void updateAfterRecalibration();

  protected:
    void registerWithParametersGuess();
    void setParameterGuess() const;
    boost::shared_ptr<SmileSection_t<T> >
    smileSection(Time optionTime, Time swapLength,
                 const Cube &sabrParametersCube) const;
    Cube sabrCalibration(const Cube &marketVolCube) const;
    void fillVolatilityCube() const;
    void createSparseSmiles() const;
    std::vector<T> spreadVolInterpolation(const Date &atmOptionDate,
                                          const Period &atmSwapTenor) const;

  private:
    mutable Cube marketVolCube_;
    mutable Cube volCubeAtmCalibrated_;
    mutable Cube sparseParameters_;
    mutable Cube denseParameters_;
    mutable std::vector<std::vector<boost::shared_ptr<SmileSection_t<T> > > >
        sparseSmiles_;
    std::vector<std::vector<Handle<Quote_t<T> > > > parametersGuessQuotes_;
    mutable Cube parametersGuess_;
    std::vector<bool> isParameterFixed_;
    bool isAtmCalibrated_;
    const boost::shared_ptr<EndCriteria_t<T>> endCriteria_;
    T maxErrorTolerance_;
    const boost::shared_ptr<OptimizationMethod_t<T>> optMethod_;
    T errorAccept_;
    const bool useMaxError_;
    const Size maxGuesses_;
    const bool backwardFlat_;
    const T cutoffStrike_;

    class PrivateObserver : public Observer {
      public:
        PrivateObserver(SwaptionVolCube1x_t<Model, T> *v) : v_(v) {}
        void update() {
            v_->setParameterGuess();
            v_->update();
        }

      private:
        SwaptionVolCube1x_t<Model, T> *v_;
    };

    boost::shared_ptr<PrivateObserver> privateObserver_;
};

//=======================================================================//
//                        SwaptionVolCube1x                              //
//=======================================================================//

template <template <class> class Model, class T>
SwaptionVolCube1x_t<Model, T>::SwaptionVolCube1x_t(
    const Handle<SwaptionVolatilityStructure_t<T>> &atmVolStructure,
    const std::vector<Period> &optionTenors,
    const std::vector<Period> &swapTenors, const std::vector<T> &strikeSpreads,
    const std::vector<std::vector<Handle<Quote_t<T> > > > &volSpreads,
    const boost::shared_ptr<SwapIndex_t<T> > &swapIndexBase,
    const boost::shared_ptr<SwapIndex_t<T> > &shortSwapIndexBase,
    bool vegaWeightedSmileFit,
    const std::vector<std::vector<Handle<Quote_t<T> > > > &parametersGuess,
    const std::vector<bool> &isParameterFixed, bool isAtmCalibrated,
    const boost::shared_ptr<EndCriteria_t<T>> &endCriteria, T maxErrorTolerance,
    const boost::shared_ptr<OptimizationMethod_t<T>> &optMethod, const T errorAccept,
    const bool useMaxError, const Size maxGuesses, const bool backwardFlat,
    const T cutoffStrike)
    : SwaptionVolatilityCube(atmVolStructure, optionTenors, swapTenors,
                             strikeSpreads, volSpreads, swapIndexBase,
                             shortSwapIndexBase, vegaWeightedSmileFit),
      parametersGuessQuotes_(parametersGuess),
      isParameterFixed_(isParameterFixed), isAtmCalibrated_(isAtmCalibrated),
      endCriteria_(endCriteria), optMethod_(optMethod),
      useMaxError_(useMaxError), maxGuesses_(maxGuesses),
      backwardFlat_(backwardFlat), cutoffStrike_(cutoffStrike) {

    if (maxErrorTolerance != Null<T>()) {
        maxErrorTolerance_ = maxErrorTolerance;
    } else {
        maxErrorTolerance_ = SWAPTIONVOLCUBE_TOL;
        if (vegaWeightedSmileFit_)
            maxErrorTolerance_ = SWAPTIONVOLCUBE_VEGAWEIGHTED_TOL;
    }
    if (errorAccept != Null<T>()) {
        errorAccept_ = errorAccept;
    } else {
        errorAccept_ = maxErrorTolerance_ / 5.0;
    }

    privateObserver_ = boost::make_shared<PrivateObserver>(this);
    this->registerWithParametersGuess();
    setParameterGuess();
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::registerWithParametersGuess() {
    for (Size i = 0; i < 4; i++)
        for (Size j = 0; j < nOptionTenors_; j++)
            for (Size k = 0; k < nSwapTenors_; k++)
                privateObserver_->registerWith(
                    parametersGuessQuotes_[j + k * nOptionTenors_][i]);
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::setParameterGuess() const {

    //! set parametersGuess_ by parametersGuessQuotes_
    parametersGuess_ = Cube(optionDates_, swapTenors_, optionTimes_,
                            swapLengths_, 4, true, backwardFlat_);
    Size i;
    for (i = 0; i < 4; i++)
        for (Size j = 0; j < nOptionTenors_; j++)
            for (Size k = 0; k < nSwapTenors_; k++) {
                parametersGuess_.setElement(
                    i, j, k,
                    parametersGuessQuotes_[j + k * nOptionTenors_][i]->value());
            }
    parametersGuess_.updateInterpolators();
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::performCalculations() const {

    SwaptionVolatilityDiscrete::performCalculations();

    //! set marketVolCube_ by volSpreads_ quotes
    marketVolCube_ =
        Cube(optionDates_, swapTenors_, optionTimes_, swapLengths_, nStrikes_);
    T atmForward;
    T atmVol, vol;
    for (Size j = 0; j < nOptionTenors_; ++j) {
        for (Size k = 0; k < nSwapTenors_; ++k) {
            atmForward = atmStrike(optionDates_[j], swapTenors_[k]);
            atmVol = atmVol_->volatility(optionDates_[j], swapTenors_[k],
                                         atmForward);
            for (Size i = 0; i < nStrikes_; ++i) {
                vol = atmVol + volSpreads_[j * nSwapTenors_ + k][i]->value();
                marketVolCube_.setElement(i, j, k, vol);
            }
        }
    }
    marketVolCube_.updateInterpolators();

    sparseParameters_ = sabrCalibration(marketVolCube_);
    // parametersGuess_ = sparseParameters_;
    sparseParameters_.updateInterpolators();
    // parametersGuess_.updateInterpolators();
    volCubeAtmCalibrated_ = marketVolCube_;

    if (isAtmCalibrated_) {
        fillVolatilityCube();
        denseParameters_ = sabrCalibration(volCubeAtmCalibrated_);
        denseParameters_.updateInterpolators();
    }
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::updateAfterRecalibration() {
    volCubeAtmCalibrated_ = marketVolCube_;
    if (isAtmCalibrated_) {
        fillVolatilityCube();
        denseParameters_ = sabrCalibration(volCubeAtmCalibrated_);
        denseParameters_.updateInterpolators();
    }
    notifyObservers();
}

template <template <class> class Model, class T>
typename SwaptionVolCube1x_t<Model, T>::Cube
SwaptionVolCube1x_t<Model, T>::sabrCalibration(
    const Cube &marketVolCube) const {

    const std::vector<Time> &optionTimes = marketVolCube.optionTimes();
    const std::vector<Time> &swapLengths = marketVolCube.swapLengths();
    const std::vector<Date> &optionDates = marketVolCube.optionDates();
    const std::vector<Period> &swapTenors = marketVolCube.swapTenors();
    Matrix_t<T> alphas(optionTimes.size(), swapLengths.size(), 0.);
    Matrix_t<T> betas(alphas);
    Matrix_t<T> nus(alphas);
    Matrix_t<T> rhos(alphas);
    Matrix_t<T> forwards(alphas);
    Matrix_t<T> errors(alphas);
    Matrix_t<T> maxErrors(alphas);
    Matrix_t<T> endCriteria(alphas);

    const std::vector<Matrix_t<T> > &tmpMarketVolCube = marketVolCube.points();

    std::vector<T> strikes(strikeSpreads_.size());
    std::vector<T> volatilities(strikeSpreads_.size());

    for (Size j = 0; j < optionTimes.size(); j++) {
        for (Size k = 0; k < swapLengths.size(); k++) {
            T atmForward = atmStrike(optionDates[j], swapTenors[k]);
            strikes.clear();
            volatilities.clear();
            for (Size i = 0; i < nStrikes_; i++) {
                T strike = atmForward + strikeSpreads_[i];
                if (strike >= cutoffStrike_) {
                    strikes.push_back(strike);
                    volatilities.push_back(tmpMarketVolCube[i][j][k]);
                }
            }

            const std::vector<T> &guess =
                parametersGuess_.operator()(optionTimes[j], swapLengths[k]);

            const boost::shared_ptr<typename Model<T>::Interpolation>
                sabrInterpolation =
                    boost::shared_ptr<typename Model<T>::Interpolation>(
                        new (typename Model<T>::Interpolation)(
                            strikes.begin(), strikes.end(),
                            volatilities.begin(), optionTimes[j], atmForward,
                            guess[0], guess[1], guess[2], guess[3],
                            isParameterFixed_[0], isParameterFixed_[1],
                            isParameterFixed_[2], isParameterFixed_[3],
                            vegaWeightedSmileFit_, endCriteria_, optMethod_,
                            errorAccept_, useMaxError_, maxGuesses_));
            sabrInterpolation->update();

            T rmsError = sabrInterpolation->rmsError();
            T maxError = sabrInterpolation->maxError();
            alphas[j][k] = sabrInterpolation->alpha();
            betas[j][k] = sabrInterpolation->beta();
            nus[j][k] = sabrInterpolation->nu();
            rhos[j][k] = sabrInterpolation->rho();
            forwards[j][k] = atmForward;
            errors[j][k] = rmsError;
            maxErrors[j][k] = maxError;
            endCriteria[j][k] = sabrInterpolation->endCriteria();

            QL_ENSURE(endCriteria[j][k] != EndCriteria_t<T>::MaxIterations,
                      "global swaptions calibration failed: "
                      "MaxIterations reached: "
                          << "\n"
                          << "option maturity = " << optionDates[j] << ", \n"
                          << "swap tenor = " << swapTenors[k] << ", \n"
                          << "error = " << io::rate(errors[j][k]) << ", \n"
                          << "max error = " << io::rate(maxErrors[j][k])
                          << ", \n"
                          << "   alpha = " << alphas[j][k] << "n"
                          << "   beta = " << betas[j][k] << "\n"
                          << "   nu = " << nus[j][k] << "\n"
                          << "   rho = " << rhos[j][k] << "\n");

            QL_ENSURE(useMaxError_ ? maxError : rmsError < maxErrorTolerance_,
                      "global swaptions calibration failed: "
                      "option tenor "
                          << optionDates[j] << ", swap tenor " << swapTenors[k]
                          << (useMaxError_ ? ": max error " : ": error")
                          << (useMaxError_ ? maxError : rmsError)
                          << "   alpha = " << alphas[j][k] << "n"
                          << "   beta = " << betas[j][k] << "\n"
                          << "   nu = " << nus[j][k] << "\n"
                          << "   rho = " << rhos[j][k] << "\n"
                          << (useMaxError_ ? ": error" : ": max error ")
                          << (useMaxError_ ? rmsError : maxError));
        }
    }
    Cube sabrParametersCube(optionDates, swapTenors, optionTimes, swapLengths,
                            8, true, backwardFlat_);
    sabrParametersCube.setLayer(0, alphas);
    sabrParametersCube.setLayer(1, betas);
    sabrParametersCube.setLayer(2, nus);
    sabrParametersCube.setLayer(3, rhos);
    sabrParametersCube.setLayer(4, forwards);
    sabrParametersCube.setLayer(5, errors);
    sabrParametersCube.setLayer(6, maxErrors);
    sabrParametersCube.setLayer(7, endCriteria);

    return sabrParametersCube;
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::sabrCalibrationSection(
    const Cube &marketVolCube, Cube &parametersCube,
    const Period &swapTenor) const {

    const std::vector<Time> &optionTimes = marketVolCube.optionTimes();
    const std::vector<Time> &swapLengths = marketVolCube.swapLengths();
    const std::vector<Date> &optionDates = marketVolCube.optionDates();
    const std::vector<Period> &swapTenors = marketVolCube.swapTenors();

    Size k = std::find(swapTenors.begin(), swapTenors.end(), swapTenor) -
             swapTenors.begin();
    QL_REQUIRE(k != swapTenors.size(), "swap tenor not found");

    std::vector<T> calibrationResult(8, 0.);
    const std::vector<Matrix_t<T> > &tmpMarketVolCube = marketVolCube.points();

    std::vector<T> strikes(strikeSpreads_.size());
    std::vector<T> volatilities(strikeSpreads_.size());

    for (Size j = 0; j < optionTimes.size(); j++) {
        T atmForward = atmStrike(optionDates[j], swapTenors[k]);
        strikes.clear();
        volatilities.clear();
        for (Size i = 0; i < nStrikes_; i++) {
            T strike = atmForward + strikeSpreads_[i];
            if (strike >= cutoffStrike_) {
                strikes.push_back(strike);
                volatilities.push_back(tmpMarketVolCube[i][j][k]);
            }
        }

        const std::vector<T> &guess =
            parametersGuess_.operator()(optionTimes[j], swapLengths[k]);

        const boost::shared_ptr<SABRInterpolation_t<T> > sabrInterpolation =
            boost::shared_ptr<SABRInterpolation_t<T> >(
                new SABRInterpolation_t<T>(
                    strikes.begin(), strikes.end(), volatilities.begin(),
                    optionTimes[j], atmForward, guess[0], guess[1], guess[2],
                    guess[3], isParameterFixed_[0], isParameterFixed_[1],
                    isParameterFixed_[2], isParameterFixed_[3],
                    vegaWeightedSmileFit_, endCriteria_, optMethod_,
                    errorAccept_, useMaxError_, maxGuesses_));

        sabrInterpolation->update();
        T interpolationError = sabrInterpolation->rmsError();
        calibrationResult[0] = sabrInterpolation->alpha();
        calibrationResult[1] = sabrInterpolation->beta();
        calibrationResult[2] = sabrInterpolation->nu();
        calibrationResult[3] = sabrInterpolation->rho();
        calibrationResult[4] = atmForward;
        calibrationResult[5] = interpolationError;
        calibrationResult[6] = sabrInterpolation->maxError();
        calibrationResult[7] = sabrInterpolation->endCriteria();

        QL_ENSURE(
            calibrationResult[7] != EndCriteria_t<T>::MaxIterations,
            "section calibration failed: "
            "option tenor "
                << optionDates[j] << ", swap tenor " << swapTenors[k]
                << ": max iteration (" << endCriteria_->maxIterations() << ")"
                << ", alpha " << calibrationResult[0] << ", beta "
                << calibrationResult[1] << ", nu " << calibrationResult[2]
                << ", rho " << calibrationResult[3] << ", max error "
                << calibrationResult[6] << ", error " << calibrationResult[5]);

        QL_ENSURE(
            useMaxError_ ? calibrationResult[6]
                         : calibrationResult[5] < maxErrorTolerance_,
            "section calibration failed: "
            "option tenor "
                << optionDates[j] << ", swap tenor " << swapTenors[k]
                << (useMaxError_ ? ": max error " : ": error ")
                << (useMaxError_ ? calibrationResult[6] : calibrationResult[5])
                << ", alpha " << calibrationResult[0] << ", beta "
                << calibrationResult[1] << ", nu " << calibrationResult[2]
                << ", rho " << calibrationResult[3]
                << (useMaxError_ ? ": error" : ": max error ")
                << (useMaxError_ ? calibrationResult[5]
                                 : calibrationResult[6]));

        parametersCube.setPoint(optionDates[j], swapTenors[k], optionTimes[j],
                                swapLengths[k], calibrationResult);
        parametersCube.updateInterpolators();
    }
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::fillVolatilityCube() const {

    const boost::shared_ptr<SwaptionVolatilityDiscrete> atmVolStructure =
        boost::dynamic_pointer_cast<SwaptionVolatilityDiscrete>(*atmVol_);

    std::vector<Time> atmOptionTimes(atmVolStructure->optionTimes());
    std::vector<Time> optionTimes(volCubeAtmCalibrated_.optionTimes());
    atmOptionTimes.insert(atmOptionTimes.end(), optionTimes.begin(),
                          optionTimes.end());
    std::sort(atmOptionTimes.begin(), atmOptionTimes.end());
    std::vector<Time>::iterator new_end =
        std::unique(atmOptionTimes.begin(), atmOptionTimes.end());
    atmOptionTimes.erase(new_end, atmOptionTimes.end());

    std::vector<Time> atmSwapLengths(atmVolStructure->swapLengths());
    std::vector<Time> swapLengths(volCubeAtmCalibrated_.swapLengths());
    atmSwapLengths.insert(atmSwapLengths.end(), swapLengths.begin(),
                          swapLengths.end());
    std::sort(atmSwapLengths.begin(), atmSwapLengths.end());
    new_end = std::unique(atmSwapLengths.begin(), atmSwapLengths.end());
    atmSwapLengths.erase(new_end, atmSwapLengths.end());

    std::vector<Date> atmOptionDates = atmVolStructure->optionDates();
    std::vector<Date> optionDates(volCubeAtmCalibrated_.optionDates());
    atmOptionDates.insert(atmOptionDates.end(), optionDates.begin(),
                          optionDates.end());
    std::sort(atmOptionDates.begin(), atmOptionDates.end());
    std::vector<Date>::iterator new_end_1 =
        std::unique(atmOptionDates.begin(), atmOptionDates.end());
    atmOptionDates.erase(new_end_1, atmOptionDates.end());

    std::vector<Period> atmSwapTenors = atmVolStructure->swapTenors();
    std::vector<Period> swapTenors(volCubeAtmCalibrated_.swapTenors());
    atmSwapTenors.insert(atmSwapTenors.end(), swapTenors.begin(),
                         swapTenors.end());
    std::sort(atmSwapTenors.begin(), atmSwapTenors.end());
    std::vector<Period>::iterator new_end_2 =
        std::unique(atmSwapTenors.begin(), atmSwapTenors.end());
    atmSwapTenors.erase(new_end_2, atmSwapTenors.end());

    createSparseSmiles();

    for (Size j = 0; j < atmOptionTimes.size(); j++) {

        for (Size k = 0; k < atmSwapLengths.size(); k++) {
            bool expandOptionTimes = !(std::binary_search(
                optionTimes.begin(), optionTimes.end(), atmOptionTimes[j]));
            bool expandSwapLengths = !(std::binary_search(
                swapLengths.begin(), swapLengths.end(), atmSwapLengths[k]));
            if (expandOptionTimes || expandSwapLengths) {
                T atmForward = atmStrike(atmOptionDates[j], atmSwapTenors[k]);
                T atmVol = atmVol_->volatility(atmOptionDates[j],
                                               atmSwapTenors[k], atmForward);
                std::vector<T> spreadVols =
                    spreadVolInterpolation(atmOptionDates[j], atmSwapTenors[k]);
                std::vector<T> volAtmCalibrated;
                volAtmCalibrated.reserve(nStrikes_);
                for (Size i = 0; i < nStrikes_; i++)
                    volAtmCalibrated.push_back(atmVol + spreadVols[i]);
                volCubeAtmCalibrated_.setPoint(
                    atmOptionDates[j], atmSwapTenors[k], atmOptionTimes[j],
                    atmSwapLengths[k], volAtmCalibrated);
            }
        }
    }
    volCubeAtmCalibrated_.updateInterpolators();
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::createSparseSmiles() const {

    std::vector<Time> optionTimes(sparseParameters_.optionTimes());
    std::vector<Time> swapLengths(sparseParameters_.swapLengths());
    sparseSmiles_.clear();

    for (Size j = 0; j < optionTimes.size(); j++) {
        std::vector<boost::shared_ptr<SmileSection_t<T> > > tmp;
        Size n = swapLengths.size();
        tmp.reserve(n);
        for (Size k = 0; k < n; ++k) {
            tmp.push_back(smileSection(optionTimes[j], swapLengths[k],
                                       sparseParameters_));
        }
        sparseSmiles_.push_back(tmp);
    }
}

template <template <class> class Model, class T>
std::vector<T> SwaptionVolCube1x_t<Model, T>::spreadVolInterpolation(
    const Date &atmOptionDate, const Period &atmSwapTenor) const {

    Time atmOptionTime = timeFromReference(atmOptionDate);
    Time atmTimeLength = swapLength(atmSwapTenor);

    std::vector<T> result;
    const std::vector<Time> &optionTimes(sparseParameters_.optionTimes());
    const std::vector<Time> &swapLengths(sparseParameters_.swapLengths());
    const std::vector<Date> &optionDates = sparseParameters_.optionDates();
    const std::vector<Period> &swapTenors = sparseParameters_.swapTenors();

    typename std::vector<T>::const_iterator optionTimesPreviousNode,
        swapLengthsPreviousNode;

    optionTimesPreviousNode =
        std::lower_bound(optionTimes.begin(), optionTimes.end(), atmOptionTime);
    Size optionTimesPreviousIndex =
        optionTimesPreviousNode - optionTimes.begin();
    if (optionTimesPreviousIndex > 0)
        optionTimesPreviousIndex--;

    swapLengthsPreviousNode =
        std::lower_bound(swapLengths.begin(), swapLengths.end(), atmTimeLength);
    Size swapLengthsPreviousIndex =
        swapLengthsPreviousNode - swapLengths.begin();
    if (swapLengthsPreviousIndex > 0)
        swapLengthsPreviousIndex--;

    std::vector<std::vector<boost::shared_ptr<SmileSection_t<T> > > > smiles;
    std::vector<boost::shared_ptr<SmileSection_t<T> > > smilesOnPreviousExpiry;
    std::vector<boost::shared_ptr<SmileSection_t<T> > > smilesOnNextExpiry;

    QL_REQUIRE(optionTimesPreviousIndex + 1 < sparseSmiles_.size(),
               "optionTimesPreviousIndex+1 >= sparseSmiles_.size()");
    QL_REQUIRE(swapLengthsPreviousIndex + 1 < sparseSmiles_[0].size(),
               "swapLengthsPreviousIndex+1 >= sparseSmiles_[0].size()");
    smilesOnPreviousExpiry.push_back(
        sparseSmiles_[optionTimesPreviousIndex][swapLengthsPreviousIndex]);
    smilesOnPreviousExpiry.push_back(
        sparseSmiles_[optionTimesPreviousIndex][swapLengthsPreviousIndex + 1]);
    smilesOnNextExpiry.push_back(
        sparseSmiles_[optionTimesPreviousIndex + 1][swapLengthsPreviousIndex]);
    smilesOnNextExpiry.push_back(
        sparseSmiles_[optionTimesPreviousIndex + 1][swapLengthsPreviousIndex +
                                                    1]);

    smiles.push_back(smilesOnPreviousExpiry);
    smiles.push_back(smilesOnNextExpiry);

    std::vector<T> optionsNodes(2);
    optionsNodes[0] = optionTimes[optionTimesPreviousIndex];
    optionsNodes[1] = optionTimes[optionTimesPreviousIndex + 1];

    std::vector<Date> optionsDateNodes(2);
    optionsDateNodes[0] = optionDates[optionTimesPreviousIndex];
    optionsDateNodes[1] = optionDates[optionTimesPreviousIndex + 1];

    std::vector<T> swapLengthsNodes(2);
    swapLengthsNodes[0] = swapLengths[swapLengthsPreviousIndex];
    swapLengthsNodes[1] = swapLengths[swapLengthsPreviousIndex + 1];

    std::vector<Period> swapTenorNodes(2);
    swapTenorNodes[0] = swapTenors[swapLengthsPreviousIndex];
    swapTenorNodes[1] = swapTenors[swapLengthsPreviousIndex + 1];

    T atmForward = atmStrike(atmOptionDate, atmSwapTenor);

    Matrix_t<T> atmForwards(2, 2, 0.0);
    Matrix_t<T> atmVols(2, 2, 0.0);
    for (Size i = 0; i < 2; i++) {
        for (Size j = 0; j < 2; j++) {
            atmForwards[i][j] =
                atmStrike(optionsDateNodes[i], swapTenorNodes[j]);
            // atmVols[i][j] = smiles[i][j]->volatility(atmForwards[i][j]);
            atmVols[i][j] = atmVol_->volatility(
                optionsDateNodes[i], swapTenorNodes[j], atmForwards[i][j]);
            /* With the old implementation the interpolated spreads on ATM
               volatilities were null even if the spreads on ATM volatilities to
               be
               interpolated were non-zero. The new implementation removes
               this behaviour, but introduces a small ERROR in the cube:
               even if no spreads are applied on any cube ATM volatility
               corresponding
               to quoted smile sections (that is ATM volatilities in sparse
               cube), the
               cube ATM volatilities corresponding to not quoted smile sections
               (that
               is ATM volatilities in dense cube) are no more exactly the quoted
               values,
               but that ones PLUS the linear interpolation of the fit errors on
               the ATM
               volatilities in sparse cube whose spreads are used in the
               calculation.
               A similar imprecision is introduced to the volatilities in dense
               cube
               whith moneyness near to 1.
               (See below how spreadVols are calculated).
               The extent of this error depends on the quality of the fit: in
               case of
               good fits it is negligibile.
            */
        }
    }

    for (Size k = 0; k < nStrikes_; k++) {
        const T strike =
            QLFCT::max(atmForward + strikeSpreads_[k], cutoffStrike_);
        const T moneyness = atmForward / strike;

        Matrix_t<T> strikes(2, 2, 0.);
        Matrix_t<T> spreadVols(2, 2, 0.);
        for (Size i = 0; i < 2; i++) {
            for (Size j = 0; j < 2; j++) {
                strikes[i][j] = atmForwards[i][j] / moneyness;
                spreadVols[i][j] =
                    smiles[i][j]->volatility(strikes[i][j]) - atmVols[i][j];
            }
        }
        Cube localInterpolator(optionsDateNodes, swapTenorNodes, optionsNodes,
                               swapLengthsNodes, 1);
        localInterpolator.setLayer(0, spreadVols);
        localInterpolator.updateInterpolators();

        result.push_back(localInterpolator(atmOptionTime, atmTimeLength)[0]);
    }
    return result;
}

template <template <class> class Model, class T>
boost::shared_ptr<SmileSection_t<T> >
SwaptionVolCube1x_t<Model, T>::smileSection(
    Time optionTime, Time swapLength, const Cube &sabrParametersCube) const {

    this->calculate();
    const std::vector<T> sabrParameters =
        sabrParametersCube(optionTime, swapLength);
    return boost::shared_ptr<SmileSection_t<T>>(
        new (typename Model<T>::SmileSection)(optionTime, sabrParameters[4],
                                              sabrParameters));
}

template <template <class> class Model, class T>
boost::shared_ptr<SmileSection_t<T> >
SwaptionVolCube1x_t<Model, T>::smileSectionImpl(Time optionTime,
                                                Time swapLength) const {
    if (isAtmCalibrated_)
        return smileSection(optionTime, swapLength, denseParameters_);
    else
        return smileSection(optionTime, swapLength, sparseParameters_);
}

template <template <class> class Model, class T>
Matrix_t<T> SwaptionVolCube1x_t<Model, T>::sparseSabrParameters() const {
    this->calculate();
    return sparseParameters_.browse();
}

template <template <class> class Model, class T>
Matrix_t<T> SwaptionVolCube1x_t<Model, T>::denseSabrParameters() const {
    this->calculate();
    return denseParameters_.browse();
}

template <template <class> class Model, class T>
Matrix_t<T> SwaptionVolCube1x_t<Model, T>::marketVolCube() const {
    this->calculate();
    return marketVolCube_.browse();
}

template <template <class> class Model, class T>
Matrix_t<T> SwaptionVolCube1x_t<Model, T>::volCubeAtmCalibrated() const {
    this->calculate();
    return volCubeAtmCalibrated_.browse();
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::recalibration(T beta,
                                                  const Period &swapTenor) {

    std::vector<T> betaVector(nOptionTenors_, beta);
    recalibration(betaVector, swapTenor);
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::recalibration(const std::vector<T> &beta,
                                                  const Period &swapTenor) {

    QL_REQUIRE(beta.size() == nOptionTenors_,
               "beta size (" << beta.size()
                             << ") must be equal to number of option tenors ("
                             << nOptionTenors_ << ")");

    const std::vector<Period> &swapTenors = marketVolCube_.swapTenors();
    Size k = std::find(swapTenors.begin(), swapTenors.end(), swapTenor) -
             swapTenors.begin();

    QL_REQUIRE(k != swapTenors.size(), "swap tenor (" << swapTenor
                                                      << ") not found");

    for (Size i = 0; i < nOptionTenors_; ++i) {
        parametersGuess_.setElement(1, i, k, beta[i]);
    }

    parametersGuess_.updateInterpolators();
    sabrCalibrationSection(marketVolCube_, sparseParameters_, swapTenor);

    volCubeAtmCalibrated_ = marketVolCube_;
    if (isAtmCalibrated_) {
        fillVolatilityCube();
        sabrCalibrationSection(volCubeAtmCalibrated_, denseParameters_,
                               swapTenor);
    }
    notifyObservers();
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::recalibration(
    const std::vector<Period> &swapLengths, const std::vector<T> &beta,
    const Period &swapTenor) {

    QL_REQUIRE(beta.size() == swapLengths.size(),
               "beta size (" << beta.size()
                             << ") must be equal to number of swap lenghts ("
                             << swapLengths.size() << ")");

    std::vector<Time> betaTimes;
    for (Size i = 0; i < beta.size(); i++)
        betaTimes.push_back(
            timeFromReference(optionDateFromTenor(swapLengths[i])));

    LinearInterpolation betaInterpolation(betaTimes.begin(), betaTimes.end(),
                                          beta.begin());

    std::vector<T> cubeBeta;
    for (Size i = 0; i < optionTimes().size(); i++) {
        T t = optionTimes()[i];
        // flat extrapolation ensures admissable values
        if (t < betaTimes.front())
            t = betaTimes.front();
        if (t > betaTimes.back())
            t = betaTimes.back();
        cubeBeta.push_back(betaInterpolation(t));
    }

    recalibration(cubeBeta, swapTenor);
}

//======================================================================//
//                      SwaptionVolCube1x::Cube                         //
//======================================================================//

template <template <class> class Model, class T>
SwaptionVolCube1x_t<Model, T>::Cube::Cube(const std::vector<Date> &optionDates,
                                          const std::vector<Period> &swapTenors,
                                          const std::vector<Time> &optionTimes,
                                          const std::vector<Time> &swapLengths,
                                          Size nLayers, bool extrapolation,
                                          bool backwardFlat)
    : optionTimes_(optionTimes), swapLengths_(swapLengths),
      optionDates_(optionDates), swapTenors_(swapTenors), nLayers_(nLayers),
      extrapolation_(extrapolation), backwardFlat_(backwardFlat) {

    QL_REQUIRE(optionTimes.size() > 1, "Cube::Cube(...): optionTimes.size()<2");
    QL_REQUIRE(swapLengths.size() > 1, "Cube::Cube(...): swapLengths.size()<2");

    QL_REQUIRE(optionTimes.size() == optionDates.size(),
               "Cube::Cube(...): optionTimes/optionDates mismatch");
    QL_REQUIRE(swapTenors.size() == swapLengths.size(),
               "Cube::Cube(...): swapTenors/swapLengths mismatch");

    std::vector<Matrix_t<T> > points(
        nLayers_, Matrix_t<T>(optionTimes_.size(), swapLengths_.size(), 0.0));
    for (Size k = 0; k < nLayers_; k++) {
        boost::shared_ptr<Interpolation2D_t<T> > interpolation;
        transposedPoints_.push_back(transpose(points[k]));
        if (k <= 4 && backwardFlat_)
            interpolation =
                boost::make_shared<BackwardflatLinearInterpolation_t<T> >(
                    optionTimes_.begin(), optionTimes_.end(),
                    swapLengths_.begin(), swapLengths_.end(),
                    transposedPoints_[k]);
        else
            interpolation = boost::make_shared<BilinearInterpolation_t<T> >(
                optionTimes_.begin(), optionTimes_.end(), swapLengths_.begin(),
                swapLengths_.end(), transposedPoints_[k]);
        interpolators_.push_back(boost::shared_ptr<Interpolation2D_t<T> >(
                                     new FlatExtrapolator2D_t<T>(interpolation)));
        interpolators_[k]->enableExtrapolation();
    }
    setPoints(points);
}

template <template <class> class Model, class T>
SwaptionVolCube1x_t<Model, T>::Cube::Cube(const Cube &o) {
    optionTimes_ = o.optionTimes_;
    swapLengths_ = o.swapLengths_;
    optionDates_ = o.optionDates_;
    swapTenors_ = o.swapTenors_;
    nLayers_ = o.nLayers_;
    extrapolation_ = o.extrapolation_;
    backwardFlat_ = o.backwardFlat_;
    transposedPoints_ = o.transposedPoints_;
    for (Size k = 0; k < nLayers_; ++k) {
        boost::shared_ptr<Interpolation2D_t<T> > interpolation;
        if (k <= 4 && backwardFlat_)
            interpolation = boost::make_shared<BackwardflatLinearInterpolation_t<T>>(
                optionTimes_.begin(), optionTimes_.end(), swapLengths_.begin(),
                swapLengths_.end(), transposedPoints_[k]);
        else
            interpolation = boost::make_shared<BilinearInterpolation_t<T>>(
                optionTimes_.begin(), optionTimes_.end(), swapLengths_.begin(),
                swapLengths_.end(), transposedPoints_[k]);
        interpolators_.push_back(boost::shared_ptr<Interpolation2D_t<T> >(
                                     new FlatExtrapolator2D_t<T>(interpolation)));
        interpolators_[k]->enableExtrapolation();
    }
    setPoints(o.points_);
}

template <template <class> class Model, class T>
typename SwaptionVolCube1x_t<Model, T>::Cube &
    SwaptionVolCube1x_t<Model, T>::Cube::
    operator=(const Cube &o) {
    optionTimes_ = o.optionTimes_;
    swapLengths_ = o.swapLengths_;
    optionDates_ = o.optionDates_;
    swapTenors_ = o.swapTenors_;
    nLayers_ = o.nLayers_;
    extrapolation_ = o.extrapolation_;
    backwardFlat_ = o.backwardFlat_;
    transposedPoints_ = o.transposedPoints_;
    for (Size k = 0; k < nLayers_; k++) {
        boost::shared_ptr<Interpolation2D_t<T> > interpolation;
        if (k <= 4 && backwardFlat_)
            interpolation =
                boost::make_shared<BackwardflatLinearInterpolation_t<T> >(
                    optionTimes_.begin(), optionTimes_.end(),
                    swapLengths_.begin(), swapLengths_.end(),
                    transposedPoints_[k]);
        else
            interpolation = boost::make_shared<BilinearInterpolation_t<T>>(
                optionTimes_.begin(), optionTimes_.end(), swapLengths_.begin(),
                swapLengths_.end(), transposedPoints_[k]);
        interpolators_.push_back(boost::shared_ptr<Interpolation2D_t<T> >(
                                     new FlatExtrapolator2D_t<T>(interpolation)));
        interpolators_[k]->enableExtrapolation();
    }
    setPoints(o.points_);
    return *this;
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::Cube::setElement(Size IndexOfLayer,
                                                     Size IndexOfRow,
                                                     Size IndexOfColumn, T x) {
    QL_REQUIRE(IndexOfLayer < nLayers_,
               "Cube::setElement: incompatible IndexOfLayer ");
    QL_REQUIRE(IndexOfRow < optionTimes_.size(),
               "Cube::setElement: incompatible IndexOfRow");
    QL_REQUIRE(IndexOfColumn < swapLengths_.size(),
               "Cube::setElement: incompatible IndexOfColumn");
    points_[IndexOfLayer][IndexOfRow][IndexOfColumn] = x;
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::Cube::setPoints(
    const std::vector<Matrix_t<T> > &x) {
    QL_REQUIRE(x.size() == nLayers_,
               "Cube::setPoints: incompatible number of layers ");
    QL_REQUIRE(x[0].rows() == optionTimes_.size(),
               "Cube::setPoints: incompatible size 1");
    QL_REQUIRE(x[0].columns() == swapLengths_.size(),
               "Cube::setPoints: incompatible size 2");

    points_ = x;
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::Cube::setLayer(Size i,
                                                   const Matrix_t<T> &x) {
    QL_REQUIRE(i < nLayers_, "Cube::setLayer: incompatible number of layer ");
    QL_REQUIRE(x.rows() == optionTimes_.size(),
               "Cube::setLayer: incompatible size 1");
    QL_REQUIRE(x.columns() == swapLengths_.size(),
               "Cube::setLayer: incompatible size 2");

    points_[i] = x;
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::Cube::setPoint(
    const Date &optionDate, const Period &swapTenor, const Time optionTime,
    const Time swapLength, const std::vector<T> &point) {
    const bool expandOptionTimes = !(std::binary_search(
        optionTimes_.begin(), optionTimes_.end(), optionTime));
    const bool expandSwapLengths = !(std::binary_search(
        swapLengths_.begin(), swapLengths_.end(), swapLength));

    typename std::vector<T>::const_iterator optionTimesPreviousNode,
        swapLengthsPreviousNode;

    optionTimesPreviousNode =
        std::lower_bound(optionTimes_.begin(), optionTimes_.end(), optionTime);
    Size optionTimesIndex = optionTimesPreviousNode - optionTimes_.begin();

    swapLengthsPreviousNode =
        std::lower_bound(swapLengths_.begin(), swapLengths_.end(), swapLength);
    Size swapLengthsIndex = swapLengthsPreviousNode - swapLengths_.begin();

    if (expandOptionTimes || expandSwapLengths)
        expandLayers(optionTimesIndex, expandOptionTimes, swapLengthsIndex,
                     expandSwapLengths);

    for (Size k = 0; k < nLayers_; ++k)
        points_[k][optionTimesIndex][swapLengthsIndex] = point[k];

    optionTimes_[optionTimesIndex] = optionTime;
    swapLengths_[swapLengthsIndex] = swapLength;
    optionDates_[optionTimesIndex] = optionDate;
    swapTenors_[swapLengthsIndex] = swapTenor;
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::Cube::expandLayers(Size i,
                                                       bool expandOptionTimes,
                                                       Size j,
                                                       bool expandSwapLengths) {
    QL_REQUIRE(i <= optionTimes_.size(),
               "Cube::expandLayers: incompatible size 1");
    QL_REQUIRE(j <= swapLengths_.size(),
               "Cube::expandLayers: incompatible size 2");

    if (expandOptionTimes) {
        optionTimes_.insert(optionTimes_.begin() + i, 0.);
        optionDates_.insert(optionDates_.begin() + i, Date());
    }
    if (expandSwapLengths) {
        swapLengths_.insert(swapLengths_.begin() + j, 0.);
        swapTenors_.insert(swapTenors_.begin() + j, Period());
    }

    std::vector<Matrix_t<T> > newPoints(
        nLayers_, Matrix_t<T>(optionTimes_.size(), swapLengths_.size(), 0.));

    for (Size k = 0; k < nLayers_; ++k) {
        for (Size u = 0; u < points_[k].rows(); ++u) {
            Size indexOfRow = u;
            if (u >= i && expandOptionTimes)
                indexOfRow = u + 1;
            for (Size v = 0; v < points_[k].columns(); ++v) {
                Size indexOfCol = v;
                if (v >= j && expandSwapLengths)
                    indexOfCol = v + 1;
                newPoints[k][indexOfRow][indexOfCol] = points_[k][u][v];
            }
        }
    }
    setPoints(newPoints);
}

template <template <class> class Model, class T>
const std::vector<Matrix_t<T> > &
SwaptionVolCube1x_t<Model, T>::Cube::points() const {
    return points_;
}

template <template <class> class Model, class T>
std::vector<T> SwaptionVolCube1x_t<Model, T>::Cube::
operator()(const Time optionTime, const Time swapLength) const {
    std::vector<T> result;
    for (Size k = 0; k < nLayers_; ++k)
        result.push_back(interpolators_[k]->operator()(optionTime, swapLength));
    return result;
}

template <template <class> class Model, class T>
const std::vector<Time> &
SwaptionVolCube1x_t<Model, T>::Cube::optionTimes() const {
    return optionTimes_;
}

template <template <class> class Model, class T>
const std::vector<Time> &
SwaptionVolCube1x_t<Model, T>::Cube::swapLengths() const {
    return swapLengths_;
}

template <template <class> class Model, class T>
void SwaptionVolCube1x_t<Model, T>::Cube::updateInterpolators() const {
    for (Size k = 0; k < nLayers_; ++k) {
        transposedPoints_[k] = transpose(points_[k]);
        boost::shared_ptr<Interpolation2D_t<T> > interpolation;
        if (k <= 4 && backwardFlat_)
            interpolation = boost::make_shared<BackwardflatLinearInterpolation>(
                optionTimes_.begin(), optionTimes_.end(), swapLengths_.begin(),
                swapLengths_.end(), transposedPoints_[k]);
        else
            interpolation = boost::make_shared<BilinearInterpolation_t<T>>(
                optionTimes_.begin(), optionTimes_.end(), swapLengths_.begin(),
                swapLengths_.end(), transposedPoints_[k]);
        interpolators_[k] = boost::shared_ptr<Interpolation2D_t<T> >(
            new FlatExtrapolator2D_t<T>(interpolation));
        interpolators_[k]->enableExtrapolation();
    }
}

template <template <class> class Model, class T>
Matrix_t<T> SwaptionVolCube1x_t<Model, T>::Cube::browse() const {
    Matrix_t<T> result(swapLengths_.size() * optionTimes_.size(), nLayers_ + 2,
                       0.0);
    for (Size i = 0; i < swapLengths_.size(); ++i) {
        for (Size j = 0; j < optionTimes_.size(); ++j) {
            result[i * optionTimes_.size() + j][0] = swapLengths_[i];
            result[i * optionTimes_.size() + j][1] = optionTimes_[j];
            for (Size k = 0; k < nLayers_; ++k)
                result[i * optionTimes_.size() + j][2 + k] = points_[k][j][i];
        }
    }
    return result;
}

//======================================================================//
//                      SwaptionVolCube1 (Sabr)                         //
//======================================================================//

template <class T> struct SwaptionVolCubeSabrModel {
    typedef SABRInterpolation_t<T> Interpolation;
    typedef SabrSmileSection_t<T> SmileSection;
};

template <class T> struct SwaptionVolCube1_t {
    typedef SwaptionVolCube1x_t<SwaptionVolCubeSabrModel, T> type;
};

typedef SwaptionVolCube1_t<Real>::type SwaptionVolCube1;

//======================================================================//
//                      SwaptionVolCube1a (NoArbSabr)                    //
//======================================================================//

// template <class T> struct SwaptionVolCubeNoArbSabrModel {
//     typedef NoArbSabrInterpolation_t<T> Interpolation;
//     typedef NoArbSabrSmileSection_t<T> SmileSection;
// };

// template <class T> struct SwaptionVolCube1a_t {
//     typedef SwaptionVolCube1x_t<SwaptionVolCubeNoArbSabrModel, T> type;
// };

// typedef SwaptionVolCube1a_t<Real>::type SwaptionVolCube1a;

} // namespace QuantLib

#endif
