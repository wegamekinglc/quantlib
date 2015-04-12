/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
 Copyright (C) 2007 Marco Bianchetti
 Copyright (C) 2007 François du Vignaud
 Copyright (C) 2007 Giorgio Facchinetti
 Copyright (C) 2006 Mario Pucci
 Copyright (C) 2006 StatPro Italia srl
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

/*! \file xabrinterpolation.hpp
    \brief generic interpolation class for sabr style underlying models
           like the Hagan 2002 expansion, Doust's no arbitrage sabr,
           Andreasen's zabr expansion for the masses and similar
*/

#ifndef ql_xabr_interpolation_hpp
#define ql_xabr_interpolation_hpp

#include <ql/utilities/null.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/math/optimization/method.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/math/optimization/projectedcostfunction.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/randomnumbers/haltonrsg.hpp>

namespace QuantLib {

namespace detail {

template <template <class> class Model, class T> class XABRCoeffHolder_t {
  public:
    XABRCoeffHolder_t(const Time t, const T &forward, std::vector<T> params,
                      std::vector<bool> paramIsFixed)
        : t_(t), forward_(forward), params_(params),
          paramIsFixed_(paramIsFixed.size(), false), weights_(std::vector<T>()),
          error_(Null<T>()), maxError_(Null<T>()),
          XABREndCriteria_(EndCriteria_t<T>::None) {
        QL_REQUIRE(t > 0.0, "expiry time must be positive: " << t
                                                             << " not allowed");
        QL_REQUIRE(params.size() == Model<T>().dimension(),
                   "wrong number of parameters (" << params.size()
                                                  << "), should be "
                                                  << Model<T>().dimension());
        QL_REQUIRE(paramIsFixed.size() == Model<T>().dimension(),
                   "wrong number of fixed parameters flags ("
                       << paramIsFixed.size() << "), should be "
                       << Model<T>().dimension());

        for (Size i = 0; i < params.size(); ++i) {
            if (params[i] != Null<T>())
                paramIsFixed_[i] = paramIsFixed[i];
        }
        Model<T>().defaultValues(params_, paramIsFixed_, forward_, t_);
        updateModelInstance();
    }
    virtual ~XABRCoeffHolder_t() {}

    void updateModelInstance() {
        // forward might have changed
        QL_REQUIRE(forward_ > 0.0,
                   "forward must be positive: " << forward_ << " not allowed");
        modelInstance_ = Model<T>().instance(t_, forward_, params_);
    }

    /*! Expiry, Forward */
    Time t_;
    const T &forward_;
    /*! Parameters */
    std::vector<T> params_;
    std::vector<bool> paramIsFixed_;
    std::vector<T> weights_;
    /*! Interpolation results */
    T error_, maxError_;
    typename EndCriteria_t<T>::Type XABREndCriteria_;
    /*! Model instance (if required) */
    boost::shared_ptr<typename Model<T>::type> modelInstance_;
};

template <class I1, class I2, template <class> class Model, class T>
class XABRInterpolationImpl_t
    : public Interpolation_t<T>::template templateImpl<I1, I2>,
      public XABRCoeffHolder_t<Model, T> {
  public:
    XABRInterpolationImpl_t(
        const I1 &xBegin, const I1 &xEnd, const I2 &yBegin, Time t,
        const T &forward, std::vector<T> params, std::vector<bool> paramIsFixed,
        bool vegaWeighted, const boost::shared_ptr<EndCriteria_t<T>> &endCriteria,
        const boost::shared_ptr<OptimizationMethod_t<T>> &optMethod,
        const T errorAccept, const bool useMaxError, const Size maxGuesses)
        : Interpolation_t<T>::template templateImpl<I1, I2>(xBegin, xEnd,
                                                            yBegin),
          XABRCoeffHolder_t<Model, T>(t, forward, params, paramIsFixed),
          endCriteria_(endCriteria), optMethod_(optMethod),
          errorAccept_(errorAccept), useMaxError_(useMaxError),
          maxGuesses_(maxGuesses), forward_(forward),
          vegaWeighted_(vegaWeighted) {
        // if no optimization method or endCriteria is provided, we provide one
        if (!optMethod_)
            optMethod_ = boost::shared_ptr<OptimizationMethod_t<T> >(
                new LevenbergMarquardt_t<T>(1e-8, 1e-8, 1e-8));
        // optMethod_ = boost::shared_ptr<OptimizationMethod_t<T>>(new
        //    Simplex(0.01));
        if (!endCriteria_) {
            endCriteria_ = boost::shared_ptr<EndCriteria_t<T> >(
                new EndCriteria_t<T>(60000, 100, 1e-8, 1e-8, 1e-8));
        }
        this->weights_ = std::vector<T>(xEnd - xBegin, 1.0 / (xEnd - xBegin));
    }

    void update() {

        this->updateModelInstance();

        // we should also check that y contains positive values only

        // we must update weights if it is vegaWeighted
        if (vegaWeighted_) {
            typename std::vector<T>::const_iterator x = this->xBegin_;
            typename std::vector<T>::const_iterator y = this->yBegin_;
            // std::vector<T>::iterator w = weights_.begin();
            this->weights_.clear();
            T weightsSum = 0.0;
            for (; x != this->xEnd_; ++x, ++y) {
                T stdDev = QLFCT::sqrt((*y) * (*y) * this->t_);
                this->weights_.push_back(
                    blackFormulaStdDevDerivative<T>(*x, forward_, stdDev));
                weightsSum += this->weights_.back();
            }
            // weight normalization
            typename std::vector<T>::iterator w = this->weights_.begin();
            for (; w != this->weights_.end(); ++w)
                *w /= weightsSum;
        }

        // there is nothing to optimize
        if (std::accumulate(this->paramIsFixed_.begin(),
                            this->paramIsFixed_.end(), true,
                            std::logical_and<bool>())) {
            this->error_ = interpolationError();
            this->maxError_ = interpolationMaxError();
            this->XABREndCriteria_ = EndCriteria_t<T>::None;
            return;
        } else {
            XABRError costFunction(this);

            Array_t<T> guess(Model<T>().dimension());
            for (Size i = 0; i < guess.size(); ++i)
                guess[i] = this->params_[i];

            Size iterations = 0;
            Size freeParameters = 0;
            T bestError = QL_MAX_REAL;
            Array_t<T> bestParameters;
            for (Size i = 0; i < Model<T>().dimension(); ++i)
                if (!this->paramIsFixed_[i])
                    ++freeParameters;
            HaltonRsg halton(freeParameters, 42);
            typename EndCriteria_t<T>::Type tmpEndCriteria;
            T tmpInterpolationError;

            do {

                if (iterations > 0) {
                    HaltonRsg::sample_type s = halton.nextSequence();
                    std::vector<T> sTmp(s.value.begin(), s.value.end());
                    Model<T>().guess(guess, this->paramIsFixed_, forward_,
                                     this->t_, sTmp);
                    for (Size i = 0; i < this->paramIsFixed_.size(); ++i)
                        if (this->paramIsFixed_[i])
                            guess[i] = this->params_[i];
                }

                Array_t<T> inversedTransformatedGuess(Model<T>().inverse(
                    guess, this->paramIsFixed_, this->params_, forward_));

                ProjectedCostFunction_t<T> constrainedXABRError(
                    costFunction, inversedTransformatedGuess,
                    this->paramIsFixed_);

                Array_t<T> projectedGuess(
                    constrainedXABRError.project(inversedTransformatedGuess));

                NoConstraint_t<T> constraint;
                Problem_t<T> problem(constrainedXABRError, constraint,
                                projectedGuess);
                tmpEndCriteria = optMethod_->minimize(problem, *endCriteria_);
                Array_t<T> projectedResult(problem.currentValue());
                Array_t<T> transfResult(
                    constrainedXABRError.include(projectedResult));

                Array_t<T> result = Model<T>().direct(
                    transfResult, this->paramIsFixed_, this->params_, forward_);
                tmpInterpolationError = useMaxError_ ? interpolationMaxError()
                                                     : interpolationError();

                if (tmpInterpolationError < bestError) {
                    bestError = tmpInterpolationError;
                    bestParameters = result;
                    this->XABREndCriteria_ = tmpEndCriteria;
                }

            } while (++iterations < maxGuesses_ &&
                     tmpInterpolationError > errorAccept_);

            for (Size i = 0; i < bestParameters.size(); ++i)
                this->params_[i] = bestParameters[i];

            this->error_ = interpolationError();
            this->maxError_ = interpolationMaxError();
        }
    }

    T value(T x) const {
        // see master
        // QL_REQUIRE(x > 0.0, "strike must be positive: " << x
        //                                                 << " not allowed");
        return this->modelInstance_->volatility(x);
    }

    T primitive(T) const { QL_FAIL("XABR primitive not implemented"); }
    T derivative(T) const { QL_FAIL("XABR derivative not implemented"); }
    T secondDerivative(T) const {
        QL_FAIL("XABR secondDerivative not implemented");
    }

    // calculate total squared weighted difference (L2 norm)
    T interpolationSquaredError() const {
        T error, totalError = 0.0;
        typename std::vector<T>::const_iterator x = this->xBegin_;
        typename std::vector<T>::const_iterator y = this->yBegin_;
        typename std::vector<T>::const_iterator w = this->weights_.begin();
        for (; x != this->xEnd_; ++x, ++y, ++w) {
            error = (value(*x) - *y);
            totalError += error * error * (*w);
        }
        return totalError;
    }

    // calculate weighted differences
    Disposable<Array_t<T> > interpolationErrors(const Array_t<T> &) const {
        Array_t<T> results(this->xEnd_ - this->xBegin_);
        typename std::vector<T>::const_iterator x = this->xBegin_;
        typename Array_t<T>::iterator r = results.begin();
        typename std::vector<T>::const_iterator y = this->yBegin_;
        typename std::vector<T>::const_iterator w = this->weights_.begin();
        for (; x != this->xEnd_; ++x, ++r, ++w, ++y) {
            *r = (value(*x) - *y) * QLFCT::sqrt(*w);
        }
        return results;
    }

    T interpolationError() const {
        Size n = this->xEnd_ - this->xBegin_;
        T squaredError = interpolationSquaredError();
        return QLFCT::sqrt(n * squaredError / (n == 1 ? 1 : (n - 1)));
    }

    T interpolationMaxError() const {
        T error, maxError = QL_MIN_REAL;
        I1 i = this->xBegin_;
        I2 j = this->yBegin_;
        for (; i != this->xEnd_; ++i, ++j) {
            error = QLFCT::abs(value(*i) - *j);
            maxError = QLFCT::max(maxError, error);
        }
        return maxError;
    }

  private:
    class XABRError : public CostFunction_t<T> {
      public:
        XABRError(XABRInterpolationImpl_t<I1, I2, Model, T> *xabr)
            : xabr_(xabr) {}

        T value(const Array_t<T> &x) const {
            const Array_t<T> y = Model<T>().direct(
                x, xabr_->paramIsFixed_, xabr_->params_, xabr_->forward_);
            for (Size i = 0; i < xabr_->params_.size(); ++i)
                xabr_->params_[i] = y[i];
            xabr_->updateModelInstance();
            return xabr_->interpolationSquaredError();
        }

        Disposable<Array_t<T> > values(const Array_t<T> &x) const {
            const Array_t<T> y = Model<T>().direct(
                x, xabr_->paramIsFixed_, xabr_->params_, xabr_->forward_);
            for (Size i = 0; i < xabr_->params_.size(); ++i)
                xabr_->params_[i] = y[i];
            xabr_->updateModelInstance();
            return xabr_->interpolationErrors(x);
        }

      private:
        XABRInterpolationImpl_t<I1, I2, Model, T> *xabr_;
    };
    boost::shared_ptr<EndCriteria_t<T>> endCriteria_;
    boost::shared_ptr<OptimizationMethod_t<T>> optMethod_;
    const T errorAccept_;
    const bool useMaxError_;
    const Size maxGuesses_;
    const T &forward_;
    bool vegaWeighted_;
    NoConstraint constraint_;
};

} // namespace detail
} // namespace QuantLib

#endif
