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

/*! \file gaussian1dmodel.hpp
    \brief basic interface for one factor interest rate models
*/

// uncomment to enable NTL support (see below for more details and references)
// #define GAUSS1D_ENABLE_NTL

#ifndef quantlib_gaussian1dmodel_hpp
#define quantlib_gaussian1dmodel_hpp

#include <ql/models/model.hpp>
#include <ql/models/parameter.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/time/date.hpp>
#include <ql/time/period.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/utilities/null.hpp>
#include <ql/patterns/lazyobject.hpp>

#ifdef GAUSS1D_ENABLE_NTL
#include <boost/math/bindings/rr.hpp>
#endif

#include <boost/math/special_functions/erf.hpp>
#include <boost/unordered_map.hpp>

namespace QuantLib {

/*! One factor interest rate model interface class
    The only methods that must be implemented by subclasses
    are the numeraire and zerobond methods for an input array
    of state variable values. The variable $y$ is understood
    to be the standardized (zero mean, unit variance) version
    of the model's original state variable $x$.

    NTL support may be enabled by defining GAUSS1D_ENABLE_NTL in this
    file. For details on NTL see
             http://www.shoup.net/ntl/

    \warning the variance of the state process conditional on
    $x(t)=x$ must be independent of the value of $x$

*/

template <class T>
class Gaussian1dModel_t : public TermStructureConsistentModel,
                          public LazyObject {
  public:
    const boost::shared_ptr<StochasticProcess1D_t<T> > stateProcess() const;

    const T numeraire(const Time t, const T y = 0.0,
                      const Handle<YieldTermStructure> &yts =
                          Handle<YieldTermStructure>()) const;

    const T zerobond(const Time T0, const Time t = 0.0, const T y = 0.0,
                     const Handle<YieldTermStructure> &yts =
                         Handle<YieldTermStructure>()) const;

    const T numeraire(const Date &referenceDate, const T y = 0.0,
                      const Handle<YieldTermStructure> &yts =
                          Handle<YieldTermStructure>()) const;

    const T zerobond(const Date &maturity,
                     const Date &referenceDate = Null<Date>(), const T y = 0.0,
                     const Handle<YieldTermStructure> &yts =
                         Handle<YieldTermStructure>()) const;

    const T zerobondOption(
        const Option::Type &type, const Date &expiry, const Date &valueDate,
        const Date &maturity, const T strike,
        const Date &referenceDate = Null<Date>(), const T y = 0.0,
        const Handle<YieldTermStructure> &yts = Handle<YieldTermStructure>(),
        const T yStdDevs = 7.0, const Size yGridPoints = 64,
        const bool extrapolatePayoff = true,
        const bool flatPayoffExtrapolation = false) const;

    const T forwardRate(const Date &fixing,
                        const Date &referenceDate = Null<Date>(),
                        const T y = 0.0,
                        boost::shared_ptr<IborIndex> iborIdx =
                            boost::shared_ptr<IborIndex>()) const;

    const T swapRate(const Date &fixing, const Period &tenor,
                     const Date &referenceDate = Null<Date>(), const T y = 0.0,
                     boost::shared_ptr<SwapIndex> swapIdx =
                         boost::shared_ptr<SwapIndex>()) const;

    const T swapAnnuity(const Date &fixing, const Period &tenor,
                        const Date &referenceDate = Null<Date>(),
                        const T y = 0.0,
                        boost::shared_ptr<SwapIndex> swapIdx =
                            boost::shared_ptr<SwapIndex>()) const;

    /*! Computes the integral
    \f[ {2\pi}^{-0.5} \int_{a}^{b} p(x) \exp{-0.5*x*x} \mathrm{d}x \f]
    with
    \f[ p(x) = ax^4+bx^3+cx^2+dx+e \f].
    */
    const static T gaussianPolynomialIntegral(const T a, const T b, const T c,
                                              const T d, const T e, const T x0,
                                              const T x1);

    /*! Computes the integral
    \f[ {2\pi}^{-0.5} \int_{a}^{b} p(x) \exp{-0.5*x*x} \mathrm{d}x \f]
    with
    \f[ p(x) = a(x-h)^4+b(x-h)^3+c(x-h)^2+d(x-h)+e \f].
    */
    const static T gaussianShiftedPolynomialIntegral(const T a, const T b,
                                                     const T c, const T d,
                                                     const T e, const T h,
                                                     const T x0, const T x1);

    /*! Generates a grid of values for the standardized state variable $y$
       at time $T$
        conditional on $y(t)=y$, covering yStdDevs standard deviations
       consisting of
        2*gridPoints+1 points */

    const Disposable<Array_t<T> > yGrid(const T yStdDevs, const int gridPoints,
                                        const T T0 = 1.0, const T t = 0,
                                        const T y = 0) const;

  private:
    // It is of great importance for performance reasons to cache underlying
    // swaps generated from indexes. In addition the indexes may only be given
    // as templates for the conventions with the tenor replaced by the actual
    // one later on.

    struct CachedSwapKey {
        const boost::shared_ptr<SwapIndex> index;
        const Date fixing;
        const Period tenor;
        const bool operator==(const CachedSwapKey &o) const {
            return index->name() == o.index->name() && fixing == o.fixing &&
                   tenor == o.tenor;
        }
    };

    struct CachedSwapKeyHasher
        : std::unary_function<CachedSwapKey, std::size_t> {
        std::size_t operator()(CachedSwapKey const &x) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, x.index->name());
            boost::hash_combine(seed, x.fixing.serialNumber());
            boost::hash_combine(seed, x.tenor.length());
            boost::hash_combine(seed, x.tenor.units());
            return seed;
        }
    };

    typedef boost::unordered_map<CachedSwapKey, boost::shared_ptr<VanillaSwap>,
                                 CachedSwapKeyHasher> CacheType;

    mutable CacheType swapCache_;

  protected:
    // we let derived classes register with the termstructure
    Gaussian1dModel_t(const Handle<YieldTermStructure> &yieldTermStructure)
        : TermStructureConsistentModel(yieldTermStructure) {
        registerWith(Settings::instance().evaluationDate());
    }

    virtual ~Gaussian1dModel_t() {}

    virtual const T
    numeraireImpl(const Time t, const T y,
                  const Handle<YieldTermStructure> &yts) const = 0;

    virtual const T
    zerobondImpl(const Time T0, const Time t, const T y,
                 const Handle<YieldTermStructure> &yts) const = 0;

    void performCalculations() const {
        evaluationDate_ = Settings::instance().evaluationDate();
        enforcesTodaysHistoricFixings_ =
            Settings::instance().enforcesTodaysHistoricFixings();
    }

    void generateArguments() {
        calculate();
        notifyObservers();
    }

    // retrieve underlying swap from cache if possible, otherwise
    // create it and store it in the cache
    boost::shared_ptr<VanillaSwap_t<T> >
    underlyingSwap(const boost::shared_ptr<SwapIndex_t<T> > &index,
                   const Date &expiry, const Period &tenor) const {

        CachedSwapKey k = {index, expiry, tenor};
        typename CacheType::iterator i = swapCache_.find(k);
        if (i == swapCache_.end()) {
            boost::shared_ptr<VanillaSwap_t<T> > underlying =
                index->clone(tenor)->underlyingSwap(expiry);
            swapCache_.insert(std::make_pair(k, underlying));
            return underlying;
        }
        return i->second;
    }

    boost::shared_ptr<StochasticProcess1D_t<T> > stateProcess_;
    mutable Date evaluationDate_;
    mutable bool enforcesTodaysHistoricFixings_;
};

typedef Gaussian1dModel_t<Real> Gaussian1dModel;

// inline definitions

template <class T>
inline const boost::shared_ptr<StochasticProcess1D_t<T> >
Gaussian1dModel_t<T>::stateProcess() const {

    QL_REQUIRE(stateProcess_ != NULL, "state process not set");
    return stateProcess_;
}

template <class T>
inline const T
Gaussian1dModel_t<T>::numeraire(const Time t, const T y,
                                const Handle<YieldTermStructure> &yts) const {

    return numeraireImpl(t, y, yts);
}

template <class T>
inline const T
Gaussian1dModel_t<T>::zerobond(const Time T0, const Time t, const T y,
                               const Handle<YieldTermStructure> &yts) const {
    return zerobondImpl(T0, t, y, yts);
}

template <class T>
inline const T
Gaussian1dModel_t<T>::numeraire(const Date &referenceDate, const T y,
                                const Handle<YieldTermStructure> &yts) const {

    return numeraire(termStructure()->timeFromReference(referenceDate), y, yts);
}

template <class T>
inline const T
Gaussian1dModel_t<T>::zerobond(const Date &maturity, const Date &referenceDate,
                               const T y,
                               const Handle<YieldTermStructure> &yts) const {

    return zerobond(termStructure()->timeFromReference(maturity),
                    referenceDate != Null<Date>()
                        ? termStructure()->timeFromReference(referenceDate)
                        : 0.0,
                    y, yts);
}

// implementation

template <class T>
const T
Gaussian1dModel_t<T>::forwardRate(const Date &fixing, const Date &referenceDate,
                                  const T y,
                                  boost::shared_ptr<IborIndex> iborIdx) const {

    QL_REQUIRE(iborIdx != NULL, "no ibor index given");

    calculate();

    if (fixing <= (evaluationDate_ + (enforcesTodaysHistoricFixings_ ? 0 : -1)))
        return iborIdx->fixing(fixing);

    Handle<YieldTermStructure> yts =
        iborIdx->forwardingTermStructure(); // might be empty, then use
                                            // model curve

    Date valueDate = iborIdx->valueDate(fixing);
    Date endDate = iborIdx->fixingCalendar().advance(
        valueDate, iborIdx->tenor(), iborIdx->businessDayConvention(),
        iborIdx->endOfMonth());
    // FIXME Here we should use the calculation date calendar ?
    T dcf = iborIdx->dayCounter().yearFraction(valueDate, endDate);

    return (zerobond(valueDate, referenceDate, y, yts) -
            zerobond(endDate, referenceDate, y, yts)) /
           (dcf * zerobond(endDate, referenceDate, y, yts));
}

template <class T>
const T
Gaussian1dModel_t<T>::swapRate(const Date &fixing, const Period &tenor,
                               const Date &referenceDate, const T y,
                               boost::shared_ptr<SwapIndex> swapIdx) const {

    QL_REQUIRE(swapIdx != NULL, "no swap index given");

    calculate();

    if (fixing <= (evaluationDate_ + (enforcesTodaysHistoricFixings_ ? 0 : -1)))
        return swapIdx->fixing(fixing);

    Handle<YieldTermStructure> ytsf =
        swapIdx->iborIndex()->forwardingTermStructure();
    Handle<YieldTermStructure> ytsd =
        swapIdx->discountingTermStructure(); // either might be empty, then
                                             // use model curve

    Schedule sched, floatSched;

    boost::shared_ptr<VanillaSwap> underlying =
        underlyingSwap(swapIdx, fixing, tenor);

    sched = underlying->fixedSchedule();

    boost::shared_ptr<OvernightIndexedSwapIndex> oisIdx =
        boost::dynamic_pointer_cast<OvernightIndexedSwapIndex>(swapIdx);
    if (oisIdx != NULL) {
        floatSched = sched;
    } else {
        floatSched = underlying->floatingSchedule();
    }

    T annuity = swapAnnuity(fixing, tenor, referenceDate, y,
                            swapIdx); // should be fine for
                                      // overnightindexed swap indices as
                                      // well
    T floatleg = 0.0;
    if (ytsf.empty() && ytsd.empty()) { // simple 100-formula can be used
                                        // only in one curve setup
        floatleg =
            (zerobond(sched.dates().front(), referenceDate, y) -
             zerobond(sched.calendar().adjust(sched.dates().back(),
                                              underlying->paymentConvention()),
                      referenceDate, y));
    } else {
        for (Size i = 1; i < floatSched.size(); i++) {
            floatleg +=
                (zerobond(floatSched[i - 1], referenceDate, y, ytsf) /
                     zerobond(floatSched[i], referenceDate, y, ytsf) -
                 1.0) *
                zerobond(floatSched.calendar().adjust(
                             floatSched[i], underlying->paymentConvention()),
                         referenceDate, y, ytsd);
        }
    }
    return floatleg / annuity;
}

template <class T>
const T
Gaussian1dModel_t<T>::swapAnnuity(const Date &fixing, const Period &tenor,
                                  const Date &referenceDate, const T y,
                                  boost::shared_ptr<SwapIndex> swapIdx) const {

    QL_REQUIRE(swapIdx != NULL, "no swap index given");

    calculate();

    Handle<YieldTermStructure> ytsd =
        swapIdx->discountingTermStructure(); // might be empty, then use
                                             // model curve

    boost::shared_ptr<VanillaSwap> underlying =
        underlyingSwap(swapIdx, fixing, tenor);

    Schedule sched = underlying->fixedSchedule();

    T annuity = 0.0;
    for (unsigned int j = 1; j < sched.size(); j++) {
        annuity += zerobond(sched.calendar().adjust(
                                sched.date(j), underlying->paymentConvention()),
                            referenceDate, y, ytsd) *
                   swapIdx->dayCounter().yearFraction(sched.date(j - 1),
                                                      sched.date(j));
    }
    return annuity;
}

template <class T>
const T Gaussian1dModel_t<T>::zerobondOption(
    const Option::Type &type, const Date &expiry, const Date &valueDate,
    const Date &maturity, const T strike, const Date &referenceDate, const T y,
    const Handle<YieldTermStructure> &yts, const T yStdDevs,
    const Size yGridPoints, const bool extrapolatePayoff,
    const bool flatPayoffExtrapolation) const {

    calculate();

    Time fixingTime = termStructure()->timeFromReference(expiry);
    Time referenceTime =
        referenceDate == Null<Date>()
            ? 0.0
            : termStructure()->timeFromReference(referenceDate);

    Array_t<T> yg = yGrid(yStdDevs, yGridPoints, fixingTime, referenceTime, y);
    Array_t<T> z = yGrid(yStdDevs, yGridPoints);

    Array_t<T> p(yg.size());

    for (Size i = 0; i < yg.size(); i++) {
        T expValDsc = zerobond(valueDate, expiry, yg[i], yts);
        T discount = zerobond(maturity, expiry, yg[i], yts) / expValDsc;
        p[i] = QLFCT::max((type == Option::Call ? 1.0 : -1.0) *
                              (discount - strike),
                          0.0) /
               numeraire(fixingTime, yg[i], yts) * expValDsc;
    }

    CubicInterpolation payoff(
        z.begin(), z.end(), p.begin(), CubicInterpolation::Spline, true,
        CubicInterpolation::Lagrange, 0.0, CubicInterpolation::Lagrange, 0.0);

    T price = 0.0;
    for (Size i = 0; i < z.size() - 1; i++) {
        price += gaussianShiftedPolynomialIntegral(
            0.0, payoff.cCoefficients()[i], payoff.bCoefficients()[i],
            payoff.aCoefficients()[i], p[i], z[i], z[i], z[i + 1]);
    }
    if (extrapolatePayoff) {
        if (flatPayoffExtrapolation) {
            price += gaussianShiftedPolynomialIntegral(
                0.0, 0.0, 0.0, 0.0, p[z.size() - 2], z[z.size() - 2],
                z[z.size() - 1], 100.0);
            price += gaussianShiftedPolynomialIntegral(0.0, 0.0, 0.0, 0.0, p[0],
                                                       z[0], -100.0, z[0]);
        } else {
            if (type == Option::Call)
                price += gaussianShiftedPolynomialIntegral(
                    0.0, payoff.cCoefficients()[z.size() - 2],
                    payoff.bCoefficients()[z.size() - 2],
                    payoff.aCoefficients()[z.size() - 2], p[z.size() - 2],
                    z[z.size() - 2], z[z.size() - 1], 100.0);
            if (type == Option::Put)
                price += gaussianShiftedPolynomialIntegral(
                    0.0, payoff.cCoefficients()[0], payoff.bCoefficients()[0],
                    payoff.aCoefficients()[0], p[0], z[0], -100.0, z[0]);
        }
    }

    return numeraire(referenceTime, y, yts) * price;
}

template <class T>
const T Gaussian1dModel_t<T>::gaussianPolynomialIntegral(const T a, const T b,
                                                         const T c, const T d,
                                                         const T e, const T y0,
                                                         const T y1) {

#ifdef GAUSS1D_ENABLE_NTL
    const boost::math::ntl::RR aa = 4.0 * a, ba = 2.0 * M_SQRT2 * b,
                               ca = 2.0 * c, da = M_SQRT2 * d;
    const boost::math::ntl::RR x0 = y0 * M_SQRT1_2, x1 = y1 * M_SQRT1_2;
    const boost::math::ntl::RR res =
        (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x1) -
         1.0 / (4.0 * M_SQRTPI) * exp(-x1 * x1) *
             (2.0 * aa * x1 * x1 * x1 + 3.0 * aa * x1 +
              2.0 * ba * (x1 * x1 + 1.0) + 2.0 * ca * x1 + 2.0 * da)) -
        (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x0) -
         1.0 / (4.0 * M_SQRTPI) * exp(-x0 * x0) *
             (2.0 * aa * x0 * x0 * x0 + 3.0 * aa * x0 +
              2.0 * ba * (x0 * x0 + 1.0) + 2.0 * ca * x0 + 2.0 * da));
    return NTL::to_double(res.value());
#else
    const T aa = 4.0 * a, ba = 2.0 * M_SQRT2 * b, ca = 2.0 * c,
            da = M_SQRT2 * d;
    const T x0 = y0 * M_SQRT1_2, x1 = y1 * M_SQRT1_2;
    return (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x1) -
            1.0 / (4.0 * M_SQRTPI) * QLFCT::exp(-x1 * x1) *
                (2.0 * aa * x1 * x1 * x1 + 3.0 * aa * x1 +
                 2.0 * ba * (x1 * x1 + 1.0) + 2.0 * ca * x1 + 2.0 * da)) -
           (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x0) -
            1.0 / (4.0 * M_SQRTPI) * QLFCT::exp(-x0 * x0) *
                (2.0 * aa * x0 * x0 * x0 + 3.0 * aa * x0 +
                 2.0 * ba * (x0 * x0 + 1.0) + 2.0 * ca * x0 + 2.0 * da));
#endif
}

template <class T>
const T Gaussian1dModel_t<T>::gaussianShiftedPolynomialIntegral(
    const T a, const T b, const T c, const T d, const T e, const T h,
    const T x0, const T x1) {
    return gaussianPolynomialIntegral(
        a, -4.0 * a * h + b, 6.0 * a * h * h - 3.0 * b * h + c,
        -4 * a * h * h * h + 3.0 * b * h * h - 2.0 * c * h + d,
        a * h * h * h * h - b * h * h * h + c * h * h - d * h + e, x0, x1);
}

template <class T>
const Disposable<Array_t<T> >
Gaussian1dModel_t<T>::yGrid(const T stdDevs, const int gridPoints, const T T0,
                            const T t, const T y) const {

    // we use that the standard deviation is independent of $x$ here !

    QL_REQUIRE(stateProcess_ != NULL, "state process not set");

    Array_t<T> result(2 * gridPoints + 1, 0.0);

    T x_t, e_0_t, e_t_T, stdDev_0_t, stdDev_t_T;
    T stdDev_0_T = stateProcess_->stdDeviation(0.0, 0.0, T0);
    T e_0_T = stateProcess_->expectation(0.0, 0.0, T0);

    if (t < QL_EPSILON) {
        // stdDev_0_t = 0.0;
        stdDev_t_T = stdDev_0_T;
        // e_0_t = 0.0;
        // x_t = 0.0;
        e_t_T = e_0_T;
    } else {
        stdDev_0_t = stateProcess_->stdDeviation(0.0, 0.0, t);
        stdDev_t_T = stateProcess_->stdDeviation(t, 0.0, T0 - t);
        e_0_t = stateProcess_->expectation(0.0, 0.0, t);
        x_t = y * stdDev_0_t + e_0_t;
        e_t_T = stateProcess_->expectation(t, x_t, T0 - t);
    }

    T h = stdDevs / ((Real)gridPoints);

    for (int j = -gridPoints; j <= gridPoints; j++) {
        result[j + gridPoints] =
            (e_t_T + stdDev_t_T * ((Real)j) * h - e_0_T) / stdDev_0_T;
    }

    return result;
}
}

#endif
