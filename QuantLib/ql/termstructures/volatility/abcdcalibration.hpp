/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2006 Cristina Duminuco
 Copyright (C) 2007 Giorgio Facchinetti
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

#ifndef quantlib_abcdcalibration_hpp
#define quantlib_abcdcalibration_hpp

#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/projectedcostfunction.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/array.hpp>
#include <ql/quote.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace QuantLib {

class OptimizationMethod;
class ParametersTransformation;

template <class T> class AbcdCalibration_t {

  private:
    class AbcdError : public CostFunction {
      public:
        AbcdError(AbcdCalibration_t<T> *abcd) : abcd_(abcd) {}

        T value(const Array &x) const {
            const Array_t<T> y = abcd_->transformation_->direct(x);
            abcd_->a_ = y[0];
            abcd_->b_ = y[1];
            abcd_->c_ = y[2];
            abcd_->d_ = y[3];
            return abcd_->error();
        }
        Disposable<Array_t<T> > values(const Array &x) const {
            const Array_t<T> y = abcd_->transformation_->direct(x);
            abcd_->a_ = y[0];
            abcd_->b_ = y[1];
            abcd_->c_ = y[2];
            abcd_->d_ = y[3];
            return abcd_->errors();
        }

      private:
        AbcdCalibration_t<T> *abcd_;
    };

    class AbcdParametersTransformation : public ParametersTransformation {
        mutable Array_t<T> y_;
        const T eps1_;

      public:
        AbcdParametersTransformation() : y_(Array_t<T>(4)), eps1_(.000000001) {}

        Array_t<T> direct(const Array &x) const {
            y_[0] = x[0] * x[0] - x[3] * x[3] + eps1_; // a + d > 0
            y_[1] = x[1];
            y_[2] = x[2] * x[2] + eps1_; // c > 0
            y_[3] = x[3] * x[3] + eps1_; // d > 0
            return y_;
        }

        Array_t<T> inverse(const Array &x) const {
            y_[0] = QLFCT::sqrt(x[0] + x[3] - eps1_);
            y_[1] = x[1];
            y_[2] = QLFCT::sqrt(x[2] - eps1_);
            y_[3] = QLFCT::sqrt(x[3] - eps1_);
            return y_;
        }
    };

  public:
    AbcdCalibration_t(){};
    AbcdCalibration_t(const std::vector<T> &t, const std::vector<T> &blackVols,
                      T aGuess = -0.06, T bGuess = 0.17, T cGuess = 0.54,
                      T dGuess = 0.17, bool aIsFixed = false,
                      bool bIsFixed = false, bool cIsFixed = false,
                      bool dIsFixed = false, bool vegaWeighted = false,
                      const boost::shared_ptr<EndCriteria> &endCriteria =
                          boost::shared_ptr<EndCriteria>(),
                      const boost::shared_ptr<OptimizationMethod> &method =
                          boost::shared_ptr<OptimizationMethod>());

    //! adjustment factors needed to match Black vols
    std::vector<T> k(const std::vector<T> &t,
                     const std::vector<T> &blackVols) const;
    void compute();
    // calibration results
    T value(T x) const;
    T error() const;
    T maxError() const;
    Disposable<Array_t<T> > errors() const;
    EndCriteria::Type endCriteria() const;

    T a() const;
    T b() const;
    T c() const;
    T d() const;

    bool aIsFixed_, bIsFixed_, cIsFixed_, dIsFixed_;
    T a_, b_, c_, d_;
    boost::shared_ptr<ParametersTransformation> transformation_;

  private:
    // optimization method used for fitting
    mutable EndCriteria::Type abcdEndCriteria_;
    boost::shared_ptr<EndCriteria> endCriteria_;
    boost::shared_ptr<OptimizationMethod> optMethod_;
    mutable std::vector<T> weights_;
    bool vegaWeighted_;
    //! Parameters
    std::vector<T> times_, blackVols_;
};

typedef AbcdCalibration_t<Real> AbcdCalibration;

// implementation

template <class T>
AbcdCalibration_t<T>::AbcdCalibration_t(
    const std::vector<T> &t, const std::vector<T> &blackVols, T a, T b, T c,
    T d, bool aIsFixed, bool bIsFixed, bool cIsFixed, bool dIsFixed,
    bool vegaWeighted, const boost::shared_ptr<EndCriteria> &endCriteria,
    const boost::shared_ptr<OptimizationMethod> &optMethod)
    : aIsFixed_(aIsFixed), bIsFixed_(bIsFixed), cIsFixed_(cIsFixed),
      dIsFixed_(dIsFixed), a_(a), b_(b), c_(c), d_(d),
      abcdEndCriteria_(EndCriteria::None), endCriteria_(endCriteria),
      optMethod_(optMethod), weights_(blackVols.size(), 1.0 / blackVols.size()),
      vegaWeighted_(vegaWeighted), times_(t), blackVols_(blackVols) {

    QL_REQUIRE(blackVols.size() == t.size(),
               "mismatch between number of times ("
                   << t.size() << ") and blackVols (" << blackVols.size()
                   << ")");

    // if no optimization method or endCriteria is provided, we provide one
    if (!optMethod_)
        optMethod_ = boost::shared_ptr<OptimizationMethod>(
            new LevenbergMarquardt(1e-8, 1e-8, 1e-8));
    // method_ = boost::shared_ptr<OptimizationMethod>(new
    //    Simplex(0.01));
    if (!endCriteria_)
        // endCriteria_ = boost::shared_ptr<EndCriteria>(new
        //    EndCriteria(60000, 100, 1e-8, 1e-8, 1e-8));
        endCriteria_ = boost::shared_ptr<EndCriteria>(
            new EndCriteria(1000, 100, 1.0e-8, 0.3e-4, 0.3e-4)); // Why 0.3e-4 ?
}

template <class T> void AbcdCalibration_t<T>::compute() {
    if (vegaWeighted_) {
        T weightsSum = 0.0;
        for (Size i = 0; i < times_.size(); i++) {
            T stdDev = QLFCT::sqrt(blackVols_[i] * blackVols_[i] * times_[i]);
            // when strike==forward, the blackFormulaStdDevDerivative becomes
            weights_[i] =
                CumulativeNormalDistribution_t<T>().derivative(.5 * stdDev);
            weightsSum += weights_[i];
        }
        // weight normalization
        for (Size i = 0; i < times_.size(); i++) {
            weights_[i] /= weightsSum;
        }
    }

    // there is nothing to optimize
    if (aIsFixed_ && bIsFixed_ && cIsFixed_ && dIsFixed_) {
        abcdEndCriteria_ = EndCriteria::None;
        // error_ = interpolationError();
        // maxError_ = interpolationMaxError();
        return;
    } else {

        AbcdError costFunction(this);
        transformation_ = boost::shared_ptr<ParametersTransformation>(
            new AbcdParametersTransformation);

        Array_t<T> guess(4);
        guess[0] = a_;
        guess[1] = b_;
        guess[2] = c_;
        guess[3] = d_;

        std::vector<bool> parameterAreFixed(4);
        parameterAreFixed[0] = aIsFixed_;
        parameterAreFixed[1] = bIsFixed_;
        parameterAreFixed[2] = cIsFixed_;
        parameterAreFixed[3] = dIsFixed_;

        Array_t<T> inversedTransformatedGuess(transformation_->inverse(guess));

        ProjectedCostFunction projectedAbcdCostFunction(
            costFunction, inversedTransformatedGuess, parameterAreFixed);

        Array_t<T> projectedGuess(
            projectedAbcdCostFunction.project(inversedTransformatedGuess));

        NoConstraint constraint;
        Problem problem(projectedAbcdCostFunction, constraint, projectedGuess);
        abcdEndCriteria_ = optMethod_->minimize(problem, *endCriteria_);
        Array_t<T> projectedResult(problem.currentValue());
        Array_t<T> transfResult(
            projectedAbcdCostFunction.include(projectedResult));

        Array_t<T> result = transformation_->direct(transfResult);
        a_ = result[0];
        b_ = result[1];
        c_ = result[2];
        d_ = result[3];

        validateAbcdParameters(a_, b_, c_, d_);
    }
}

template <class T> T AbcdCalibration_t<T>::a() const { return a_; }

template <class T> T AbcdCalibration_t<T>::b() const { return b_; }

template <class T> T AbcdCalibration_t<T>::c() const { return c_; }

template <class T> T AbcdCalibration_t<T>::d() const { return d_; }

template <class T> T AbcdCalibration_t<T>::value(T x) const {
    return abcdBlackVolatility(x, a_, b_, c_, d_);
}

template <class T>
std::vector<T> AbcdCalibration_t<T>::k(const std::vector<T> &t,
                                       const std::vector<T> &blackVols) const {
    QL_REQUIRE(blackVols.size() == t.size(),
               "mismatch between number of times ("
                   << t.size() << ") and blackVols (" << blackVols.size()
                   << ")");
    std::vector<T> k(t.size());
    for (Size i = 0; i < t.size(); i++) {
        k[i] = blackVols[i] / value(t[i]);
    }
    return k;
}

template <class T> T AbcdCalibration_t<T>::error() const {
    Size n = times_.size();
    T error, squaredError = 0.0;
    for (Size i = 0; i < times_.size(); i++) {
        error = (value(times_[i]) - blackVols_[i]);
        squaredError += error * error * (weights_[i]);
    }
    return QLFCT::sqrt(n * squaredError / (n - 1));
}

template <class T> T AbcdCalibration_t<T>::maxError() const {
    T error, maxError = QL_MIN_REAL;
    for (Size i = 0; i < times_.size(); i++) {
        error = QLFCT::abs(value(times_[i]) - blackVols_[i]);
        maxError = QLFCT::max(maxError, error);
    }
    return maxError;
}

// calculate weighted differences
template <class T>
Disposable<Array_t<T> > AbcdCalibration_t<T>::errors() const {
    Array_t<T> results(times_.size());
    for (Size i = 0; i < times_.size(); i++) {
        results[i] =
            (value(times_[i]) - blackVols_[i]) * QLFCT::sqrt(weights_[i]);
    }
    return results;
}

template <class T> EndCriteria::Type AbcdCalibration_t<T>::endCriteria() const {
    return abcdEndCriteria_;
}
}

#endif
