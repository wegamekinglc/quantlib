/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
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

/*! \file calibrationhelper.hpp
    \brief Calibration helper class
*/

#ifndef quantlib_interest_rate_modelling_calibration_helper_h
#define quantlib_interest_rate_modelling_calibration_helper_h

#include <ql/quote.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <ql/models/calibrationhelper.hpp>
#include <ql/math/solvers1d/brent.hpp>

#include <list>

namespace QuantLib {

class PricingEngine;

//! liquid market instrument used during calibration
template <class T> class CalibrationHelper_t : public LazyObject {
  public:
    enum CalibrationErrorType {
        RelativePriceError,
        PriceError,
        ImpliedVolError
    };

    CalibrationHelper_t(
        const Handle<Quote_t<T> > &volatility,
        const Handle<YieldTermStructure_t<T> > &termStructure,
        CalibrationErrorType calibrationErrorType = RelativePriceError)
        : volatility_(volatility), termStructure_(termStructure),
          calibrationErrorType_(calibrationErrorType) {
        registerWith(volatility_);
        registerWith(termStructure_);
    }

    void performCalculations() const {
        marketValue_ = blackPrice(volatility_->value());
    }

    //! returns the volatility Handle
    Handle<Quote_t<T> > volatility() { return volatility_; }

    //! returns the actual price of the instrument (from volatility)
    T marketValue() const {
        calculate();
        return marketValue_;
    }

    //! returns the price of the instrument according to the model
    virtual T modelValue() const = 0;

    //! returns the error resulting from the model valuation
    virtual T calibrationError();

    virtual void addTimesTo(std::list<Time> &times) const = 0;

    //! Black volatility implied by the model
    T impliedVolatility(T targetValue, T accuracy, Size maxEvaluations,
                        T minVol, T maxVol) const;

    //! Black price given a volatility
    virtual T blackPrice(T volatility) const = 0;

    void setPricingEngine(const boost::shared_ptr<PricingEngine> &engine) {
        engine_ = engine;
    }

  protected:
    mutable T marketValue_;
    Handle<Quote_t<T> > volatility_;
    Handle<YieldTermStructure_t<T> > termStructure_;
    boost::shared_ptr<PricingEngine> engine_;

  private:
    class ImpliedVolatilityHelper;
    const CalibrationErrorType calibrationErrorType_;
};

typedef CalibrationHelper_t<Real> CalibrationHelper;

// implementation

template <class T> class CalibrationHelper_t<T>::ImpliedVolatilityHelper {
  public:
    ImpliedVolatilityHelper(const CalibrationHelper_t<T> &helper, T value)
        : helper_(helper), value_(value) {}

    T operator()(T x) const { return value_ - helper_.blackPrice(x); }

  private:
    const CalibrationHelper_t &helper_;
    T value_;
};

template <class T>
T CalibrationHelper_t<T>::impliedVolatility(T targetValue, T accuracy,
                                            Size maxEvaluations, T minVol,
                                            T maxVol) const {

    ImpliedVolatilityHelper f(*this, targetValue);
    Brent solver;
    solver.setMaxEvaluations(maxEvaluations);
    return solver.solve(f, accuracy, volatility_->value(), minVol, maxVol);
}

template <class T> T CalibrationHelper_t<T>::calibrationError() {
    T error;

    switch (calibrationErrorType_) {
    case RelativePriceError:
        error = QLFCT::abs(marketValue() - modelValue()) / marketValue();
        break;
    case PriceError:
        error = marketValue() - modelValue();
        break;
    case ImpliedVolError: {
        const T lowerPrice = blackPrice(0.001);
        const T upperPrice = blackPrice(10);
        const T modelPrice = modelValue();

        T implied;
        if (modelPrice <= lowerPrice)
            implied = 0.001;
        else if (modelPrice >= upperPrice)
            implied = 10.0;
        else
            implied =
                this->impliedVolatility(modelPrice, 1e-12, 5000, 0.001, 10);
        error = implied - volatility_->value();
    } break;
    default:
        QL_FAIL("unknown Calibration Error Type");
    }

    return error;
}
}

#endif
