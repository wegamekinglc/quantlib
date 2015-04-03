/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2006 François du Vignaud
 Copyright (C) 2006, 2008 Ferdinando Ametrano
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

/*! \file swaptionvolmatrix.hpp
    \brief Swaption at-the-money volatility matrix
*/

#ifndef quantlib_swaption_volatility_matrix_hpp
#define quantlib_swaption_volatility_matrix_hpp

#include <ql/termstructures/volatility/swaption/swaptionvoldiscrete.hpp>
#include <ql/math/interpolations/interpolation2d.hpp>
#include <ql/math/matrix.hpp>
#include <ql/quote.hpp>
#include <ql/termstructures/volatility/flatsmilesection.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/math/interpolations/bilinearinterpolation.hpp>

#include <boost/noncopyable.hpp>
#include <vector>

namespace QuantLib {

//! At-the-money swaption-volatility matrix
/*! This class provides the at-the-money volatility for a given
    swaption by interpolating a volatility matrix whose elements
    are the market volatilities of a set of swaption with given
    option date and swapLength.

    The volatility matrix <tt>M</tt> must be defined so that:
    - the number of rows equals the number of option dates;
    - the number of columns equals the number of swap tenors;
    - <tt>M[i][j]</tt> contains the volatility corresponding
      to the <tt>i</tt>-th option and <tt>j</tt>-th tenor.
*/
template <class T>
class SwaptionVolatilityMatrix_t : public SwaptionVolatilityDiscrete,
                                   private boost::noncopyable {
  public:
    //! floating reference date, floating market data
    SwaptionVolatilityMatrix_t(
        const Calendar &calendar, BusinessDayConvention bdc,
        const std::vector<Period> &optionTenors,
        const std::vector<Period> &swapTenors,
        const std::vector<std::vector<Handle<Quote_t<T> > > > &vols,
        const DayCounter &dayCounter);
    //! fixed reference date, floating market data
    SwaptionVolatilityMatrix_t(
        const Date &referenceDate, const Calendar &calendar,
        BusinessDayConvention bdc, const std::vector<Period> &optionTenors,
        const std::vector<Period> &swapTenors,
        const std::vector<std::vector<Handle<Quote_t<T> > > > &vols,
        const DayCounter &dayCounter);
    //! floating reference date, fixed market data
    SwaptionVolatilityMatrix_t(const Calendar &calendar,
                               BusinessDayConvention bdc,
                               const std::vector<Period> &optionTenors,
                               const std::vector<Period> &swapTenors,
                               const Matrix_t<T> &volatilities,
                               const DayCounter &dayCounter);
    //! fixed reference date, fixed market data
    SwaptionVolatilityMatrix_t(const Date &referenceDate,
                               const Calendar &calendar,
                               BusinessDayConvention bdc,
                               const std::vector<Period> &optionTenors,
                               const std::vector<Period> &swapTenors,
                               const Matrix_t<T> &volatilities,
                               const DayCounter &dayCounter);
    // fixed reference date and fixed market data, option dates
    SwaptionVolatilityMatrix_t(const Date &referenceDate,
                               const std::vector<Date> &optionDates,
                               const std::vector<Period> &swapTenors,
                               const Matrix_t<T> &volatilities,
                               const DayCounter &dayCounter);
    //! \name LazyObject interface
    //@{
    void performCalculations() const;
    //@}
    //! \name TermStructure interface
    //@{
    Date maxDate() const;
    //@}
    //! \name VolatilityTermStructure interface
    //@{
    T minStrike() const;
    T maxStrike() const;
    //@}
    //! \name SwaptionVolatilityStructure interface
    //@{
    const Period &maxSwapTenor() const;
    //@}
    //! \name Other inspectors
    //@{
    //! returns the lower indexes of surrounding volatility matrix corners
    std::pair<Size, Size> locate(const Date &optionDate,
                                 const Period &swapTenor) const {
        return locate(timeFromReference(optionDate), swapLength(swapTenor));
    }
    //! returns the lower indexes of surrounding volatility matrix corners
    std::pair<Size, Size> locate(Time optionTime, Time swapLength) const {
        return std::make_pair(interpolation_.locateY(optionTime),
                              interpolation_.locateX(swapLength));
    }
    //@}
  protected:
    // defining the following method would break CMS test suite
    // to be further investigated
    // boost::shared_ptr<SmileSection_t<T>> smileSectionImpl(const Date&,
    //                                                 const Period&) const;
    boost::shared_ptr<SmileSection_t<T> > smileSectionImpl(Time, Time) const;
    T volatilityImpl(Time optionTime, Time swapLength, T strike) const;

  private:
    void checkInputs(Size volRows, Size volsColumns) const;
    void registerWithMarketData();
    std::vector<std::vector<Handle<Quote_t<T> > > > volHandles_;
    mutable Matrix_t<T> volatilities_;
    Interpolation2D_t<T> interpolation_;
};

// inline definitions

template <class T> inline Date SwaptionVolatilityMatrix_t<T>::maxDate() const {
    return optionDates_.back();
}

template <class T> inline T SwaptionVolatilityMatrix_t<T>::minStrike() const {
    return QL_MIN_REAL;
}

template <class T> inline T SwaptionVolatilityMatrix_t<T>::maxStrike() const {
    return QL_MAX_REAL;
}

template <class T>
inline const Period &SwaptionVolatilityMatrix_t<T>::maxSwapTenor() const {
    return swapTenors_.back();
}

template <class T>
inline T SwaptionVolatilityMatrix_t<T>::volatilityImpl(Time optionTime,
                                                       Time swapLength,
                                                       T) const {
    calculate();
    return interpolation_(swapLength, optionTime, true);
}

typedef SwaptionVolatilityMatrix_t<Real> SwaptionVolatilityMatrix;

// implementation

// floating reference date, floating market data
template <class T>
SwaptionVolatilityMatrix_t<T>::SwaptionVolatilityMatrix_t(
    const Calendar &cal, BusinessDayConvention bdc,
    const std::vector<Period> &optionT, const std::vector<Period> &swapT,
    const std::vector<std::vector<Handle<Quote_t<T> > > > &vols,
    const DayCounter &dc)
    : SwaptionVolatilityDiscrete(optionT, swapT, 0, cal, bdc, dc),
      volHandles_(vols), volatilities_(vols.size(), vols.front().size()) {
    checkInputs(volatilities_.rows(), volatilities_.columns());
    registerWithMarketData();
    interpolation_ = BilinearInterpolation_t<T>(
        swapLengths_.begin(), swapLengths_.end(), optionTimes_.begin(),
        optionTimes_.end(), volatilities_);
}

// fixed reference date, floating market data
template <class T>
SwaptionVolatilityMatrix_t<T>::SwaptionVolatilityMatrix_t(
    const Date &refDate, const Calendar &cal, BusinessDayConvention bdc,
    const std::vector<Period> &optionT, const std::vector<Period> &swapT,
    const std::vector<std::vector<Handle<Quote_t<T> > > > &vols,
    const DayCounter &dc)
    : SwaptionVolatilityDiscrete(optionT, swapT, refDate, cal, bdc, dc),
      volHandles_(vols), volatilities_(vols.size(), vols.front().size()) {
    checkInputs(volatilities_.rows(), volatilities_.columns());
    registerWithMarketData();
    interpolation_ = BilinearInterpolation_t<T>(
        swapLengths_.begin(), swapLengths_.end(), optionTimes_.begin(),
        optionTimes_.end(), volatilities_);
}

// floating reference date, fixed market data
template <class T>
SwaptionVolatilityMatrix_t<T>::SwaptionVolatilityMatrix_t(
    const Calendar &cal, BusinessDayConvention bdc,
    const std::vector<Period> &optionT, const std::vector<Period> &swapT,
    const Matrix_t<T> &vols, const DayCounter &dc)
    : SwaptionVolatilityDiscrete(optionT, swapT, 0, cal, bdc, dc),
      volHandles_(vols.rows()), volatilities_(vols.rows(), vols.columns()) {

    checkInputs(vols.rows(), vols.columns());

    // fill dummy handles to allow generic handle-based
    // computations later on
    for (Size i = 0; i < vols.rows(); ++i) {
        volHandles_[i].resize(vols.columns());
        for (Size j = 0; j < vols.columns(); ++j)
            volHandles_[i][j] =
                Handle<Quote_t<T> >(boost::shared_ptr<Quote_t<T> >(
                    new SimpleQuote_t<T>(vols[i][j])));
    }
    interpolation_ = BilinearInterpolation_t<T>(
        swapLengths_.begin(), swapLengths_.end(), optionTimes_.begin(),
        optionTimes_.end(), volatilities_);
}

// fixed reference date, fixed market data
template <class T>
SwaptionVolatilityMatrix_t<T>::SwaptionVolatilityMatrix_t(
    const Date &refDate, const Calendar &cal, BusinessDayConvention bdc,
    const std::vector<Period> &optionT, const std::vector<Period> &swapT,
    const Matrix_t<T> &vols, const DayCounter &dc)
    : SwaptionVolatilityDiscrete(optionT, swapT, refDate, cal, bdc, dc),
      volHandles_(vols.rows()), volatilities_(vols.rows(), vols.columns()) {

    checkInputs(vols.rows(), vols.columns());

    // fill dummy handles to allow generic handle-based
    // computations later on
    for (Size i = 0; i < vols.rows(); ++i) {
        volHandles_[i].resize(vols.columns());
        for (Size j = 0; j < vols.columns(); ++j)
            volHandles_[i][j] =
                Handle<Quote_t<T> >(boost::shared_ptr<Quote_t<T> >(
                    new SimpleQuote_t<T>(vols[i][j])));
    }
    interpolation_ = BilinearInterpolation_t<T>(
        swapLengths_.begin(), swapLengths_.end(), optionTimes_.begin(),
        optionTimes_.end(), volatilities_);
}

// fixed reference date and fixed market data, option dates
template <class T>
SwaptionVolatilityMatrix_t<T>::SwaptionVolatilityMatrix_t(
    const Date &today, const std::vector<Date> &optionDates,
    const std::vector<Period> &swapT, const Matrix_t<T> &vols,
    const DayCounter &dc)
    : SwaptionVolatilityDiscrete(optionDates, swapT, today, Calendar(),
                                 Following, dc),
      volHandles_(vols.rows()), volatilities_(vols.rows(), vols.columns()) {

    checkInputs(vols.rows(), vols.columns());

    // fill dummy handles to allow generic handle-based
    // computations later on
    for (Size i = 0; i < vols.rows(); ++i) {
        volHandles_[i].resize(vols.columns());
        for (Size j = 0; j < vols.columns(); ++j)
            volHandles_[i][j] =
                Handle<Quote_t<T> >(boost::shared_ptr<Quote_t<T> >(
                    new SimpleQuote_t<T>(vols[i][j])));
    }
    interpolation_ = BilinearInterpolation_t<T>(
        swapLengths_.begin(), swapLengths_.end(), optionTimes_.begin(),
        optionTimes_.end(), volatilities_);
}

template <class T>
void SwaptionVolatilityMatrix_t<T>::checkInputs(Size volRows,
                                                Size volsColumns) const {
    QL_REQUIRE(nOptionTenors_ == volRows,
               "mismatch between number of option dates ("
                   << nOptionTenors_ << ") and number of rows (" << volRows
                   << ") in the vol matrix");
    QL_REQUIRE(nSwapTenors_ == volsColumns,
               "mismatch between number of swap tenors ("
                   << nSwapTenors_ << ") and number of rows (" << volsColumns
                   << ") in the vol matrix");
}

template <class T>
void SwaptionVolatilityMatrix_t<T>::registerWithMarketData() {
    for (Size i = 0; i < volHandles_.size(); ++i)
        for (Size j = 0; j < volHandles_.front().size(); ++j)
            registerWith(volHandles_[i][j]);
}

template <class T>
void SwaptionVolatilityMatrix_t<T>::performCalculations() const {

    SwaptionVolatilityDiscrete::performCalculations();

    // we might use iterators here...
    for (Size i = 0; i < volatilities_.rows(); ++i)
        for (Size j = 0; j < volatilities_.columns(); ++j)
            volatilities_[i][j] = volHandles_[i][j]->value();
}

// boost::shared_ptr<SmileSection_t<T>>
// SwaptionVolatilityMatrix_t<T>::smileSectionImpl(const Date& d,
//                                           const Period& swapTenor) const {
//    Time optionTime = timeFromReference(d);
//    Time swapLength = convertSwapTenor(swapTenor);
//    // dummy strike
//    T atmVol = volatilityImpl(optionTime, swapLength, 0.05);
//    return boost::shared_ptr<SmileSection_t<T>>(new
//        FlatSmileSection_t<T>(d, atmVol, dayCounter(), referenceDate()));
//}

template <class T>
boost::shared_ptr<SmileSection_t<T> >
SwaptionVolatilityMatrix_t<T>::smileSectionImpl(Time optionTime,
                                                Time swapLength) const {
    // dummy strike
    T atmVol = volatilityImpl(optionTime, swapLength, 0.05);
    return boost::shared_ptr<SmileSection_t<T> >(
        new FlatSmileSection_t<T>(optionTime, atmVol, dayCounter()));
}
}

#endif
