/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Banca Profilo S.p.A.
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

/*! \file forwardmeasureprocess.hpp
    \brief forward-measure stochastic processes
*/

#ifndef quantlib_forward_measure_processes_hpp
#define quantlib_forward_measure_processes_hpp

#include <ql/stochasticprocess.hpp>

namespace QuantLib {

//! forward-measure stochastic process
/*! stochastic process whose dynamics are expressed in the forward
    measure.

    \ingroup processes
*/
template <class T>
class ForwardMeasureProcess_t : public StochasticProcess_t<T> {
  public:
    virtual void setForwardMeasureTime(Time);
    Time getForwardMeasureTime() const;

  protected:
    ForwardMeasureProcess_t() {}
    ForwardMeasureProcess_t(Time T0) : T0_(T0) {}
    ForwardMeasureProcess_t(const boost::shared_ptr<
        typename StochasticProcess_t<T>::discretization> &);
    Time T0_;
};

//! forward-measure 1-D stochastic process
/*! 1-D stochastic process whose dynamics are expressed in the
    forward measure.

    \ingroup processes
*/
template <class T>
class ForwardMeasureProcess1D_t : public StochasticProcess1D_t<T> {
  public:
    virtual void setForwardMeasureTime(Time);
    Time getForwardMeasureTime() const;

  protected:
    ForwardMeasureProcess1D_t() {}
    ForwardMeasureProcess1D_t(Time T0) : T0_(T0) {}
    ForwardMeasureProcess1D_t(const boost::shared_ptr<
        typename StochasticProcess_t<T>::discretization> &);
    Time T0_;
};

typedef ForwardMeasureProcess_t<Real> ForwardMeasureProcess;

// implementation

template <class T>
ForwardMeasureProcess_t<T>::ForwardMeasureProcess_t(const boost::shared_ptr<
    typename StochasticProcess_t<T>::discretization> &disc)
    : StochasticProcess(disc) {}

template <class T>
void ForwardMeasureProcess_t<T>::setForwardMeasureTime(Time T0) {
    T0_ = T0;
    this->notifyObservers();
}

template <class T>
Time ForwardMeasureProcess_t<T>::getForwardMeasureTime() const {
    return T0_;
}

// 1-D specialization

template <class T>
ForwardMeasureProcess1D_t<T>::ForwardMeasureProcess1D_t(const boost::shared_ptr<
    typename StochasticProcess_t<T>::discretization> &disc)
    : StochasticProcess1D(disc) {}

template <class T>
void ForwardMeasureProcess1D_t<T>::setForwardMeasureTime(Time T0) {
    T0_ = T0;
    this->notifyObservers();
}

template <class T>
Time ForwardMeasureProcess1D_t<T>::getForwardMeasureTime() const {
    return T0_;
}

typedef ForwardMeasureProcess1D_t<Real> ForwardMeasureProcess1D;

} // namespace QuantLib

#endif
