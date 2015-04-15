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

/*! \file gsrprocess.hpp
    \brief GSR model process with piecewise volatilities and mean reversions,
           the dynamic is expressed in some T-forward measure.
           You may provide a single value for the mean reversion, then
           it is assumed to be constant. For many grid points (like 20 and
           above) evaluation may get slow. A caching is therefore provided.
           By that the results become inconsistent as soon as the parameters
           change. In that case flushCache() must be called. To ensure correct
           calibration this is done in the generateArguments() of the GSR model
*/

#ifndef quantlib_gsr_process_hpp
#define quantlib_gsr_process_hpp

#include <ql/processes/forwardmeasureprocess.hpp>
#include <ql/math/comparison.hpp>
#include <map>

namespace QuantLib {

//! GSR stochastic process
/*! \ingroup processes */
template <class T> class GsrProcess_t : public ForwardMeasureProcess1D_t<T> {
  public:
    GsrProcess_t(const Array_t<Time> &times, const Array_t<T> &vols,
                 const Array_t<T> &reversions, const Time T0 = 60.0);
    //! \name StochasticProcess1D interface
    //@{
    T x0() const;
    T drift(Time t, T x) const;
    T diffusion(Time t, T) const;
    T expectation(Time t0, T x0, Time dt) const;
    T stdDeviation(Time t0, T x0, Time dt) const;
    T variance(Time t0, T, Time dt) const;

    T sigma(Time t) const;
    T reversion(Time t) const;
    T y(Time t) const;
    T G(Time t, Time T0, T x) const;
    //@}
    void setForwardMeasureTime(Time t) {
        flushCache();
        ForwardMeasureProcess1D_t<T>::setForwardMeasureTime(t);
    }
    void flushCache() const;

  protected:
    const Array_t<Time> &times_;
    const Array_t<T> &vols_;
    const Array_t<T> &reversions_;

  private:
    T expectationp1(Time t0, T x0,
                    Time dt) const; // expectation can be split into a x0
                                    // dependent term (p1) and an
                                    // independent term (p2)
    T expectationp2(Time t0, Time dt) const;
    const int lowerIndex(Time t) const;
    const int upperIndex(Time t) const;
    const Time time2(Size index) const;
    const Time cappedTime(Size index, Time cap = Null<Time>()) const;
    const Time flooredTime(Size index, Time floor = Null<Time>()) const;
    const T vol(Size index) const;
    const T rev(Size index) const;
    const bool revZero(Size index) const;
    mutable std::map<std::pair<T, T>, T> cache1_, cache2_, cache3_, cache5_;
    mutable std::map<T, T> cache4_;
    mutable std::vector<bool> revZero_;
};

typedef GsrProcess_t<Real> GsrProcess;

// implementation

template <class T>
GsrProcess_t<T>::GsrProcess_t(const Array_t<Time> &times,
                              const Array_t<T> &vols,
                              const Array_t<T> &reversions, const Time T0)
    : ForwardMeasureProcess1D_t<T>(T0), times_(times), vols_(vols),
      reversions_(reversions), revZero_(reversions.size(), false) {
    QL_REQUIRE(times.size() == vols.size() - 1,
               "number of volatilities ("
                   << vols.size() << ") compared to number of times ("
                   << times_.size() << " must be bigger by one");
    QL_REQUIRE(times.size() == reversions.size() - 1 || reversions.size() == 1,
               "number of reversions ("
                   << vols.size() << ") compared to number of times ("
                   << times_.size() << " must be bigger by one, or exactly "
                                       "1 reversion must be given");
    for (int i = 0; i < ((int)times.size()) - 1; i++)
        QL_REQUIRE(times[i] < times[i + 1], "times must be increasing ("
                                                << times[i] << "@" << i << " , "
                                                << times[i + 1] << "@" << i + 1
                                                << ")");
    flushCache();
}

template <class T> T GsrProcess_t<T>::x0() const { return 0.0; }

template <class T> T GsrProcess_t<T>::drift(Time t, T x) const {
    QL_REQUIRE(t <= this->getForwardMeasureTime(),
               "t (" << t << ") must not be greater than forward measure time ("
                     << this->getForwardMeasureTime() << ")");
    return y(t) -
           G(t, this->getForwardMeasureTime(), x) * vol(lowerIndex(t)) *
               vol(lowerIndex(t)) -
           rev(lowerIndex(t)) * x;
}

template <class T> T GsrProcess_t<T>::diffusion(Time t, T) const {
    QL_REQUIRE(t <= this->getForwardMeasureTime(),
               "t (" << t << ") must not be greater than forward measure time ("
                     << this->getForwardMeasureTime() << ")");
    return vol(lowerIndex(t));
}

template <class T> T GsrProcess_t<T>::expectation(Time w, T xw, Time dt) const {

    Time t = w + dt;
    QL_REQUIRE(t <= this->getForwardMeasureTime(),
               "t (" << t << ") must not be greater than forward measure time ("
                     << this->getForwardMeasureTime() << ")");

    return expectationp1(w, xw, dt) + expectationp2(w, dt);
}

template <class T> void GsrProcess_t<T>::flushCache() const {
    // this method must be called if parameters change (see the note
    // in the header), so we can ensure here that the zero reversion
    // flag is kept consistent, too
    for (int i = 0; i < (int)reversions_.size(); i++)
        // if (close(reversions_[i], 0.0))
        if (QLFCT::abs(reversions_[i]) < 1E-4)
            revZero_[i] = true;
        else
            revZero_[i] = false;
    cache1_.clear();
    cache2_.clear();
    cache3_.clear();
    cache4_.clear();
    cache5_.clear();
}

template <class T>
T GsrProcess_t<T>::expectationp1(Time w, T xw, Time dt) const {
    Time t = w + dt;
    std::pair<T, T> key;
    key = std::make_pair(w, t);
    typename std::map<std::pair<T, T>, T>::const_iterator k = cache1_.find(key);
    if (k != cache1_.end())
        return xw * (k->second);
    // A(w,t)x(w)
    T res2 = 1.0;
    for (int i = lowerIndex(w); i <= upperIndex(t) - 1; i++) {
        res2 *=
            QLFCT::exp(-rev(i) * (cappedTime(i + 1, t) - flooredTime(i, w)));
    }
    cache1_.insert(std::make_pair(key, res2));
    return res2 * xw;
}

template <class T> T GsrProcess_t<T>::expectationp2(Time w, Time dt) const {

    Time t = w + dt;

    std::pair<T, T> key;
    key = std::make_pair(w, t);
    typename std::map<std::pair<T, T>, T>::const_iterator k = cache2_.find(key);
    if (k != cache2_.end())
        return k->second;

    Time T0 = this->getForwardMeasureTime();

    T res = 0.0;

    // \int A(s,t)y(s)
    for (int k = lowerIndex(w); k <= upperIndex(t) - 1; k++) {
        // l<k
        for (int l = 0; l <= k - 1; l++) {
            T res2 = 1.0;
            // alpha_l
            res2 *= revZero(l)
                        ? vol(l) * vol(l) * (time2(l + 1) - time2(l))
                        : vol(l) * vol(l) / (2.0 * rev(l)) *
                              (1.0 - QLFCT::exp(-2.0 * rev(l) *
                                                (time2(l + 1) - time2(l))));
            // zeta_i (i>k)
            for (int i = k + 1; i <= upperIndex(t) - 1; i++)
                res2 *= QLFCT::exp(-rev(i) * (cappedTime(i + 1, t) - time2(i)));
            // beta_j (j<k)
            for (int j = l + 1; j <= k - 1; j++)
                res2 *= QLFCT::exp(-2.0 * rev(j) * (time2(j + 1) - time2(j)));
            // zeta_k beta_k
            res2 *=
                revZero(k)
                    ? 2.0 * time2(k) - flooredTime(k, w) -
                          cappedTime(k + 1, t) -
                          2.0 * (time2(k) - cappedTime(k + 1, t))
                    : (QLFCT::exp(rev(k) * (2.0 * time2(k) - flooredTime(k, w) -
                                            cappedTime(k + 1, t))) -
                       QLFCT::exp(2.0 * rev(k) *
                                  (time2(k) - cappedTime(k + 1, t)))) /
                          rev(k);
            // add to sum
            res += res2;
        }
        // l=k
        T res2 = 1.0;
        // alpha_k zeta_k
        res2 *=
            revZero(k)
                ? vol(k) * vol(k) / 4.0 *
                      (4.0 * QLFCT::pow(cappedTime(k + 1, t) - time2(k), 2.0) -
                       (QLFCT::pow(flooredTime(k, w) - 2.0 * time2(k) +
                                       cappedTime(k + 1, t),
                                   2.0) +
                        QLFCT::pow(cappedTime(k + 1, t) - flooredTime(k, w),
                                   2.0)))
                : vol(k) * vol(k) / (2.0 * rev(k) * rev(k)) *
                      (QLFCT::exp(-2.0 * rev(k) *
                                  (cappedTime(k + 1, t) - time2(k))) +
                       1.0 - (QLFCT::exp(-rev(k) *
                                         (flooredTime(k, w) - 2.0 * time2(k) +
                                          cappedTime(k + 1, t))) +
                              QLFCT::exp(-rev(k) * (cappedTime(k + 1, t) -
                                                    flooredTime(k, w)))));
        // zeta_i (i>k)
        for (int i = k + 1; i <= upperIndex(t) - 1; i++)
            res2 *= QLFCT::exp(-rev(i) * (cappedTime(i + 1, t) - time2(i)));
        // no beta_j in this case ...
        res += res2;
    }

    // int -A(s,t) \sigma^2 G(s,T)
    for (int k = lowerIndex(w); k <= upperIndex(t) - 1; k++) {
        T res2 = 0.0;
        // l>k
        for (int l = k + 1; l <= upperIndex(T0) - 1; l++) {
            T res3 = 1.0;
            // eta_l
            res3 *= revZero(l)
                        ? cappedTime(l + 1, T0) - time2(l)
                        : (1.0 - QLFCT::exp(-rev(l) * (cappedTime(l + 1, T0) -
                                                       time2(l)))) /
                              rev(l);
            // zeta_i (i>k)
            for (int i = k + 1; i <= upperIndex(t) - 1; i++)
                res3 *= QLFCT::exp(-rev(i) * (cappedTime(i + 1, t) - time2(i)));
            // gamma_j (j>k)
            for (int j = k + 1; j <= l - 1; j++)
                res3 *= QLFCT::exp(-rev(j) * (time2(j + 1) - time2(j)));
            // zeta_k gamma_k
            res3 *= revZero(k)
                        ? (cappedTime(k + 1, t) - time2(k + 1) -
                           (2.0 * flooredTime(k, w) - cappedTime(k + 1, t) -
                            time2(k + 1))) /
                              2.0
                        : (QLFCT::exp(rev(k) *
                                      (cappedTime(k + 1, t) - time2(k + 1))) -
                           QLFCT::exp(rev(k) *
                                      (2.0 * flooredTime(k, w) -
                                       cappedTime(k + 1, t) - time2(k + 1)))) /
                              (2.0 * rev(k));
            // add to sum
            res2 += res3;
        }
        // l=k
        T res3 = 1.0;
        // eta_k zeta_k
        res3 *=
            revZero(k)
                ? (-QLFCT::pow(cappedTime(k + 1, t) - cappedTime(k + 1, T0),
                               2.0) -
                   2.0 * QLFCT::pow(cappedTime(k + 1, t) - flooredTime(k, w),
                                    2.0) +
                   QLFCT::pow(2.0 * flooredTime(k, w) - cappedTime(k + 1, T0) -
                                  cappedTime(k + 1, t),
                              2.0)) /
                      4.0
                : (2.0 - QLFCT::exp(rev(k) * (cappedTime(k + 1, t) -
                                              cappedTime(k + 1, T0))) -
                   (2.0 * QLFCT::exp(-rev(k) * (cappedTime(k + 1, t) -
                                                flooredTime(k, w))) -
                    QLFCT::exp(rev(k) * (2.0 * flooredTime(k, w) -
                                         cappedTime(k + 1, T0) -
                                         cappedTime(k + 1, t))))) /
                      (2.0 * rev(k) * rev(k));
        // zeta_i (i>k)
        for (int i = k + 1; i <= upperIndex(t) - 1; i++)
            res3 *= QLFCT::exp(-rev(i) * (cappedTime(i + 1, t) - time2(i)));
        // no gamma_j in this case ...
        res2 += res3;
        // add to main accumulator
        res += -vol(k) * vol(k) * res2;
    }

    cache2_.insert(std::make_pair(key, res));

    return res;
}

template <class T>
T GsrProcess_t<T>::stdDeviation(Time t0, T x0, Time dt) const {
    return QLFCT::sqrt(variance(t0, x0, dt));
}

template <class T> T GsrProcess_t<T>::variance(Time w, T, Time dt) const {

    Time t = w + dt;
    QL_REQUIRE(t <= this->getForwardMeasureTime(),
               "t (" << t << ") must not be greater than forward measure time ("
                     << this->getForwardMeasureTime() << ")");

    std::pair<T, T> key;
    key = std::make_pair(w, t);
    typename std::map<std::pair<T, T>, T>::const_iterator k = cache3_.find(key);
    if (k != cache3_.end())
        return k->second;

    T res = 0.0;
    for (int k = lowerIndex(w); k <= upperIndex(t) - 1; k++) {
        T res2 = vol(k) * vol(k);
        // zeta_k^2
        res2 *=
            revZero(k)
                ? -(flooredTime(k, w) - cappedTime(k + 1, t))
                : (1.0 - QLFCT::exp(2.0 * rev(k) * (flooredTime(k, w) -
                                                    cappedTime(k + 1, t)))) /
                      (2.0 * rev(k));
        // zeta_i (i>k)
        for (int i = k + 1; i <= upperIndex(t) - 1; i++) {
            res2 *=
                QLFCT::exp(-2.0 * rev(i) * (cappedTime(i + 1, t) - time2(i)));
        }
        res += res2;
    }

    cache3_.insert(std::make_pair(key, res));
    return res;
}

template <class T> T GsrProcess_t<T>::sigma(Time t) const {
    return vol(lowerIndex(t));
}

template <class T> T GsrProcess_t<T>::reversion(Time t) const {
    return rev(lowerIndex(t));
}

template <class T> T GsrProcess_t<T>::y(Time t) const {

    QL_REQUIRE(t >= 0.0 && t <= this->getForwardMeasureTime(),
               "y(t) should be called with t (" << t << ") in Range [0,"
                                                << this->getForwardMeasureTime()
                                                << "].");

    T key;
    key = t;
    typename std::map<T, T>::const_iterator k = cache4_.find(key);
    if (k != cache4_.end())
        return k->second;

    T res = 0.0;
    for (int i = 0; i <= upperIndex(t) - 1; i++) {
        T res2 = 1.0;
        for (int j = i + 1; j <= upperIndex(t) - 1; j++) {
            res2 *=
                QLFCT::exp(-2.0 * rev(j) * (cappedTime(j + 1, t) - time2(j)));
        }
        res2 *= revZero(i)
                    ? vol(i) * vol(i) * (cappedTime(i + 1, t) - time2(i))
                    : (vol(i) * vol(i) / (2.0 * rev(i)) *
                       (1.0 - QLFCT::exp(-2.0 * rev(i) *
                                         (cappedTime(i + 1, t) - time2(i)))));
        res += res2;
    }

    cache4_.insert(std::make_pair(key, res));
    return res;
}

template <class T> T GsrProcess_t<T>::G(Time t, Time w, T) const {

    QL_REQUIRE(w >= t, "G(t,w) should be called with w ("
                           << w << ") not lesser than t (" << t << ")");
    QL_REQUIRE(t >= 0.0 && w <= this->getForwardMeasureTime(),
               "G(t,w) should be called with (t,w)=("
                   << t << "," << w << ") in Range [0,"
                   << this->getForwardMeasureTime() << "].");

    std::pair<T, T> key;
    key = std::make_pair(w, t);
    typename std::map<std::pair<T, T>, T>::const_iterator k = cache5_.find(key);
    if (k != cache5_.end())
        return k->second;

    T res = 0.0;
    for (int i = lowerIndex(t); i <= upperIndex(w) - 1; i++) {
        T res2 = 1.0;
        for (int j = lowerIndex(t); j <= i - 1; j++) {
            res2 *= QLFCT::exp(-rev(j) * (time2(j + 1) - flooredTime(j, t)));
        }
        res2 *= revZero(i) ? cappedTime(i + 1, w) - flooredTime(i, t)
                           : (1.0 - QLFCT::exp(-rev(i) * (cappedTime(i + 1, w) -
                                                          flooredTime(i, t)))) /
                                 rev(i);
        res += res2;
    }

    cache5_.insert(std::make_pair(key, res));
    return res;
}

template <class T> const int GsrProcess_t<T>::lowerIndex(Time t) const {
    return (const int)(std::upper_bound(times_.begin(), times_.end(), t) -
                       times_.begin());
}

template <class T> const int GsrProcess_t<T>::upperIndex(Time t) const {
    if (t < QL_EPSILON)
        return 0;
    return (const int)(std::upper_bound(times_.begin(), times_.end(),
                                        t - QL_EPSILON) -
                       times_.begin() + 1);
}

template <class T>
const Time GsrProcess_t<T>::cappedTime(Size index, Time cap) const {
    return cap != Null<Time>() ? std::min(cap, time2(index)) : time2(index);
}

template <class T>
const Time GsrProcess_t<T>::flooredTime(Size index, Time floor) const {
    return floor != Null<Time>() ? std::max(floor, time2(index)) : time2(index);
}

template <class T> const Time GsrProcess_t<T>::time2(Size index) const {
    if (index == 0)
        return 0.0;
    if (index > times_.size())
        return this
            ->getForwardMeasureTime(); // FIXME how to ensure that forward
                                       // measure time is geq all times
                                       // given
    return times_[index - 1];
}

template <class T> const T GsrProcess_t<T>::vol(Size index) const {
    if (index >= vols_.size())
        return vols_.back();
    return vols_[index];
}

template <class T> const T GsrProcess_t<T>::rev(Size index) const {
    if (index >= reversions_.size())
        return reversions_.back();
    return reversions_[index];
}

template <class T> const bool GsrProcess_t<T>::revZero(Size index) const {
    if (index >= revZero_.size())
        return revZero_.back();
    return revZero_[index];
}
}

#endif
