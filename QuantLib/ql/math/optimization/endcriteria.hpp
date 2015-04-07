/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2007 Marco Bianchetti
 Copyright (C) 2001, 2002, 2003 Nicolas Di Césaré
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

/*! \file endcriteria.hpp
    \brief Optimization criteria class
*/

#ifndef quantlib_optimization_criteria_hpp
#define quantlib_optimization_criteria_hpp

#include <ql/utilities/null.hpp>
#include <ql/errors.hpp>

#include <iosfwd>

namespace QuantLib {

//! Criteria to end optimization process:
/*! - maximum number of iterations AND minimum number of iterations around
   stationary point
    - x (independent variable) stationary point
    - y=f(x) (dependent variable) stationary point
    - stationary gradient
*/
template <class T> class EndCriteria_t {
  public:
    enum Type {
        None,
        MaxIterations,
        StationaryPoint,
        StationaryFunctionValue,
        StationaryFunctionAccuracy,
        ZeroGradientNorm,
        Unknown
    };

    //! Initialization constructor
    EndCriteria_t(Size maxIterations, Size maxStationaryStateIterations,
                  T rootEpsilon, T functionEpsilon, T gradientNormEpsilon);

    // Inspectors
    Size maxIterations() const;
    Size maxStationaryStateIterations() const;
    T rootEpsilon() const;
    T functionEpsilon() const;
    T gradientNormEpsilon() const;

    /*! Test if the number of iterations is not too big
        and if a minimum point is not reached */
    bool operator()(const Size iteration, Size &statState,
                    const bool positiveOptimization, const T fold,
                    const T normgold, const T fnew, const T normgnew,
                    typename EndCriteria_t<T>::Type &ecType) const;

    /*! Test if the number of iteration is below MaxIterations */
    bool checkMaxIterations(const Size iteration,
                            typename EndCriteria_t<T>::Type &ecType) const;
    /*! Test if the root variation is below rootEpsilon */
    bool checkStationaryPoint(const T xOld, const T xNew,
                              Size &statStateIterations,
                              typename EndCriteria_t<T>::Type &ecType) const;
    /*! Test if the function variation is below functionEpsilon */
    bool
    checkStationaryFunctionValue(const T fxOld, const T fxNew,
                                 Size &statStateIterations,
                                 typename EndCriteria_t<T>::Type &ecType) const;
    /*! Test if the function value is below functionEpsilon */
    bool checkStationaryFunctionAccuracy(
        const T f, const bool positiveOptimization,
        typename EndCriteria_t<T>::Type &ecType) const;
    /*! Test if the gradient norm variation is below gradientNormEpsilon */
    // bool checkZerGradientNormValue(const T gNormOld,
    //                               const T gNormNew,
    //                               typename EndCriteria_t<T>::Type& ecType)
    //                               const;
    /*! Test if the gradient norm value is below gradientNormEpsilon */
    bool checkZeroGradientNorm(const T gNorm,
                               typename EndCriteria_t<T>::Type &ecType) const;

  protected:
    //! Maximum number of iterations
    Size maxIterations_;
    //! Maximun number of iterations in stationary state
    mutable Size maxStationaryStateIterations_;
    //! root, function and gradient epsilons
    T rootEpsilon_, functionEpsilon_, gradientNormEpsilon_;
};

template <class T>
std::ostream &operator<<(std::ostream &out,
                         typename EndCriteria_t<T>::Type ecType);

typedef EndCriteria_t<Real> EndCriteria;

// implementation

template <class T>
EndCriteria_t<T>::EndCriteria_t(Size maxIterations,
                                Size maxStationaryStateIterations,
                                T rootEpsilon, T functionEpsilon,
                                T gradientNormEpsilon)
    : maxIterations_(maxIterations),
      maxStationaryStateIterations_(maxStationaryStateIterations),
      rootEpsilon_(rootEpsilon), functionEpsilon_(functionEpsilon),
      gradientNormEpsilon_(gradientNormEpsilon) {

    if (maxStationaryStateIterations_ == Null<Size>())
        maxStationaryStateIterations_ = QLFCT::min<double>(
            static_cast<Size>(maxIterations / 2), static_cast<Size>(100));
    QL_REQUIRE(maxStationaryStateIterations_ > 1,
               "maxStationaryStateIterations_ ("
                   << maxStationaryStateIterations_
                   << ") must be greater than one");
    QL_REQUIRE(maxStationaryStateIterations_ < maxIterations_,
               "maxStationaryStateIterations_ ("
                   << maxStationaryStateIterations_
                   << ") must be less than maxIterations_ (" << maxIterations_
                   << ")");
    if (gradientNormEpsilon_ == Null<T>())
        gradientNormEpsilon_ = functionEpsilon_;
}

template <class T>
bool EndCriteria_t<T>::checkMaxIterations(
    const Size iteration, typename EndCriteria_t<T>::Type &ecType) const {
    if (iteration < maxIterations_)
        return false;
    ecType = MaxIterations;
    return true;
}

template <class T>
bool EndCriteria_t<T>::checkStationaryPoint(
    const T xOld, const T xNew, Size &statStateIterations,
    typename EndCriteria_t<T>::Type &ecType) const {
    if (QLFCT::abs(xNew - xOld) >= rootEpsilon_) {
        statStateIterations = 0;
        return false;
    }
    ++statStateIterations;
    if (statStateIterations <= maxStationaryStateIterations_)
        return false;
    ecType = StationaryPoint;
    return true;
}

template <class T>
bool EndCriteria_t<T>::checkStationaryFunctionValue(
    const T fxOld, const T fxNew, Size &statStateIterations,
    typename EndCriteria_t<T>::Type &ecType) const {
    if (QLFCT::abs(fxNew - fxOld) >= functionEpsilon_) {
        statStateIterations = 0;
        return false;
    }
    ++statStateIterations;
    if (statStateIterations <= maxStationaryStateIterations_)
        return false;
    ecType = StationaryFunctionValue;
    return true;
}

template <class T>
bool EndCriteria_t<T>::checkStationaryFunctionAccuracy(
    const T f, const bool positiveOptimization,
    typename EndCriteria_t<T>::Type &ecType) const {
    if (!positiveOptimization)
        return false;
    if (f >= functionEpsilon_)
        return false;
    ecType = StationaryFunctionAccuracy;
    return true;
}

// bool EndCriteria_t<T>::checkZerGradientNormValue(
//                                        const T gNormOld,
//                                        const T gNormNew,
//                                        typename EndCriteria_t<T>::Type&
//                                        ecType) const {
//    if (QLFCT::abs(gNormNew-gNormOld) >= gradientNormEpsilon_)
//        return false;
//    ecType = StationaryGradient;
//    return true;
//}

template <class T>
bool EndCriteria_t<T>::checkZeroGradientNorm(
    const T gradientNorm, typename EndCriteria_t<T>::Type &ecType) const {
    if (gradientNorm >= gradientNormEpsilon_)
        return false;
    ecType = ZeroGradientNorm;
    return true;
}

template <class T>
bool EndCriteria_t<T>::
operator()(const Size iteration, Size &statStateIterations,
           const bool positiveOptimization, const T fold,
           const T, // normgold,
           const T fnew, const T normgnew,
           typename EndCriteria_t<T>::Type &ecType) const {
    return checkMaxIterations(iteration, ecType) ||
           checkStationaryFunctionValue(fold, fnew, statStateIterations,
                                        ecType) ||
           checkStationaryFunctionAccuracy(fnew, positiveOptimization,
                                           ecType) ||
           checkZeroGradientNorm(normgnew, ecType);
}

// Inspectors
template <class T> Size EndCriteria_t<T>::maxIterations() const {
    return maxIterations_;
}

template <class T> Size EndCriteria_t<T>::maxStationaryStateIterations() const {
    return maxStationaryStateIterations_;
}

template <class T> T EndCriteria_t<T>::rootEpsilon() const {
    return rootEpsilon_;
}

template <class T> T EndCriteria_t<T>::functionEpsilon() const {
    return functionEpsilon_;
}

template <class T> T EndCriteria_t<T>::gradientNormEpsilon() const {
    return gradientNormEpsilon_;
}

template <class T>
std::ostream &operator<<(std::ostream &out,
                         typename EndCriteria_t<T>::Type ec) {
    switch (ec) {
    case QuantLib::EndCriteria_t<T>::None:
        return out << "None";
    case QuantLib::EndCriteria_t<T>::MaxIterations:
        return out << "MaxIterations";
    case QuantLib::EndCriteria_t<T>::StationaryPoint:
        return out << "StationaryPoint";
    case QuantLib::EndCriteria_t<T>::StationaryFunctionValue:
        return out << "StationaryFunctionValue";
    case QuantLib::EndCriteria_t<T>::StationaryFunctionAccuracy:
        return out << "StationaryFunctionAccuracy";
    case QuantLib::EndCriteria_t<T>::ZeroGradientNorm:
        return out << "ZeroGradientNorm";
    case QuantLib::EndCriteria_t<T>::Unknown:
        return out << "Unknown";
    default:
        QL_FAIL("unknown EndCriteria::Type (" << Integer(ec) << ")");
    }
}
}

#endif
