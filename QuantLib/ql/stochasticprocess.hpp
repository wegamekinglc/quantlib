/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2004, 2005 StatPro Italia srl
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

/*! \file stochasticprocess.hpp
    \brief stochastic processes
*/

#ifndef quantlib_stochastic_process_hpp
#define quantlib_stochastic_process_hpp

#include <ql/time/date.hpp>
#include <ql/patterns/observable.hpp>
#include <ql/math/matrix.hpp>

namespace QuantLib {

//! multi-dimensional stochastic process class.
/*! This class describes a stochastic process governed by
    \f[
    d\mathrm{x}_t = \mu(t, x_t)\mathrm{d}t
                  + \sigma(t, \mathrm{x}_t) \cdot d\mathrm{W}_t.
    \f]
*/
template <class T>
class StochasticProcess_t : public Observer, public Observable {
  public:
    //! discretization of a stochastic process over a given time interval
    class discretization {
      public:
        virtual ~discretization() {}
        virtual Disposable<Array_t<T> > drift(const StochasticProcess_t<T> &,
                                              Time t0, const Array_t<T> &x0,
                                              Time dt) const = 0;
        virtual Disposable<Matrix_t<T> >
        diffusion(const StochasticProcess_t<T> &, Time t0, const Array_t<T> &x0,
                  Time dt) const = 0;
        virtual Disposable<Matrix_t<T> >
        covariance(const StochasticProcess_t<T> &, Time t0,
                   const Array_t<T> &x0, Time dt) const = 0;
    };
    virtual ~StochasticProcess_t() {}
    //! \name Stochastic process interface
    //@{
    //! returns the number of dimensions of the stochastic process
    virtual Size size() const = 0;
    //! returns the number of independent factors of the process
    virtual Size factors() const;
    //! returns the initial values of the state variables
    virtual Disposable<Array_t<T> > initialValues() const = 0;
    /*! \brief returns the drift part of the equation, i.e.,
               \f$ \mu(t, \mathrm{x}_t) \f$
    */
    virtual Disposable<Array_t<T> > drift(Time t,
                                          const Array_t<T> &x) const = 0;
    /*! \brief returns the diffusion part of the equation, i.e.
               \f$ \sigma(t, \mathrm{x}_t) \f$
    */
    virtual Disposable<Matrix_t<T> > diffusion(Time t,
                                               const Array_t<T> &x) const = 0;
    /*! returns the expectation
        \f$ E(\mathrm{x}_{t_0 + \Delta t}
            | \mathrm{x}_{t_0} = \mathrm{x}_0) \f$
        of the process after a time interval \f$ \Delta t \f$
        according to the given discretization. This method can be
        overridden in derived classes which want to hard-code a
        particular discretization.
    */
    virtual Disposable<Array_t<T> > expectation(Time t0, const Array_t<T> &x0,
                                                Time dt) const;
    /*! returns the standard deviation
        \f$ S(\mathrm{x}_{t_0 + \Delta t}
            | \mathrm{x}_{t_0} = \mathrm{x}_0) \f$
        of the process after a time interval \f$ \Delta t \f$
        according to the given discretization. This method can be
        overridden in derived classes which want to hard-code a
        particular discretization.
    */
    virtual Disposable<Matrix_t<T> > stdDeviation(Time t0, const Array_t<T> &x0,
                                                  Time dt) const;
    /*! returns the covariance
        \f$ V(\mathrm{x}_{t_0 + \Delta t}
            | \mathrm{x}_{t_0} = \mathrm{x}_0) \f$
        of the process after a time interval \f$ \Delta t \f$
        according to the given discretization. This method can be
        overridden in derived classes which want to hard-code a
        particular discretization.
    */
    virtual Disposable<Matrix_t<T> > covariance(Time t0, const Array_t<T> &x0,
                                                Time dt) const;
    /*! returns the asset value after a time interval \f$ \Delta t
        \f$ according to the given discretization. By default, it
        returns
        \f[
        E(\mathrm{x}_0,t_0,\Delta t) +
        S(\mathrm{x}_0,t_0,\Delta t) \cdot \Delta \mathrm{w}
        \f]
        where \f$ E \f$ is the expectation and \f$ S \f$ the
        standard deviation.
    */
    virtual Disposable<Array_t<T> > evolve(Time t0, const Array_t<T> &x0,
                                           Time dt, const Array_t<T> &dw) const;
    /*! applies a change to the asset value. By default, it
        returns \f$ \mathrm{x} + \Delta \mathrm{x} \f$.
    */
    virtual Disposable<Array_t<T> > apply(const Array_t<T> &x0,
                                          const Array_t<T> &dx) const;
    //@}

    //! \name utilities
    //@{
    /*! returns the time value corresponding to the given date
        in the reference system of the stochastic process.

        \note As a number of processes might not need this
              functionality, a default implementation is given
              which raises an exception.
    */
    virtual Time time(const Date &) const;
    //@}

    //! \name Observer interface
    //@{
    void update();
    //@}
  protected:
    StochasticProcess_t();
    StochasticProcess_t(const boost::shared_ptr<discretization> &);
    boost::shared_ptr<discretization> discretization_;
};

//! 1-dimensional stochastic process
/*! This class describes a stochastic process governed by
    \f[
        dx_t = \mu(t, x_t)dt + \sigma(t, x_t)dW_t.
    \f]
*/
template <class T> class StochasticProcess1D_t : public StochasticProcess_t<T> {
  public:
    //! discretization of a 1-D stochastic process
    class discretization {
      public:
        virtual ~discretization() {}
        virtual T drift(const StochasticProcess1D_t<T> &, Time t0, T x0,
                        Time dt) const = 0;
        virtual T diffusion(const StochasticProcess1D_t<T> &, Time t0, T x0,
                            Time dt) const = 0;
        virtual T variance(const StochasticProcess1D_t<T> &, Time t0, T x0,
                           Time dt) const = 0;
    };
    //! \name 1-D stochastic process interface
    //@{
    //! returns the initial value of the state variable
    virtual T x0() const = 0;
    //! returns the drift part of the equation, i.e. \f$ \mu(t, x_t) \f$
    virtual T drift(Time t, T x) const = 0;
    /*! \brief returns the diffusion part of the equation, i.e.
        \f$ \sigma(t, x_t) \f$
    */
    virtual T diffusion(Time t, T x) const = 0;
    /*! returns the expectation
        \f$ E(x_{t_0 + \Delta t} | x_{t_0} = x_0) \f$
        of the process after a time interval \f$ \Delta t \f$
        according to the given discretization. This method can be
        overridden in derived classes which want to hard-code a
        particular discretization.
    */
    virtual T expectation(Time t0, T x0, Time dt) const;
    /*! returns the standard deviation
        \f$ S(x_{t_0 + \Delta t} | x_{t_0} = x_0) \f$
        of the process after a time interval \f$ \Delta t \f$
        according to the given discretization. This method can be
        overridden in derived classes which want to hard-code a
        particular discretization.
    */
    virtual T stdDeviation(Time t0, T x0, Time dt) const;
    /*! returns the variance
        \f$ V(x_{t_0 + \Delta t} | x_{t_0} = x_0) \f$
        of the process after a time interval \f$ \Delta t \f$
        according to the given discretization. This method can be
        overridden in derived classes which want to hard-code a
        particular discretization.
    */
    virtual T variance(Time t0, T x0, Time dt) const;
    /*! returns the asset value after a time interval \f$ \Delta t
        \f$ according to the given discretization. By default, it
        returns
        \f[
        E(x_0,t_0,\Delta t) + S(x_0,t_0,\Delta t) \cdot \Delta w
        \f]
        where \f$ E \f$ is the expectation and \f$ S \f$ the
        standard deviation.
    */
    virtual T evolve(Time t0, T x0, Time dt, T dw) const;
    /*! applies a change to the asset value. By default, it
        returns \f$ x + \Delta x \f$.
    */
    virtual T apply(Real x0, T dx) const;
    //@}
  protected:
    StochasticProcess1D_t();
    StochasticProcess1D_t(const boost::shared_ptr<discretization> &);
    boost::shared_ptr<discretization> discretization_;

  private:
    // StochasticProcess interface implementation
    Size size() const;
    Disposable<Array_t<T> > initialValues() const;
    Disposable<Array_t<T> > drift(Time t, const Array_t<T> &x) const;
    Disposable<Matrix_t<T> > diffusion(Time t, const Array_t<T> &x) const;
    Disposable<Array_t<T> > expectation(Time t0, const Array_t<T> &x0,
                                        Time dt) const;
    Disposable<Matrix_t<T> > stdDeviation(Time t0, const Array_t<T> &x0,
                                          Time dt) const;
    Disposable<Matrix_t<T> > covariance(Time t0, const Array_t<T> &x0,
                                        Time dt) const;
    Disposable<Array_t<T> > evolve(Time t0, const Array_t<T> &x0, Time dt,
                                   const Array_t<T> &dw) const;
    Disposable<Array_t<T> > apply(const Array_t<T> &x0,
                                  const Array_t<T> &dx) const;
};

typedef StochasticProcess_t<Real> StochasticProcess;
typedef StochasticProcess1D_t<Real> StochasticProcess1D;

// inline definitions

template <class T> inline Size StochasticProcess1D_t<T>::size() const {
    return 1;
}

template <class T>
inline Disposable<Array_t<T> > StochasticProcess1D_t<T>::initialValues() const {
    Array_t<T> a(1, x0());
    return a;
}

template <class T>
inline Disposable<Array_t<T> >
StochasticProcess1D_t<T>::drift(Time t, const Array_t<T> &x) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(x.size() == 1, "1-D array required");
#endif
    Array_t<T> a(1, drift(t, x[0]));
    return a;
}

template <class T>
inline Disposable<Matrix_t<T> >
StochasticProcess1D_t<T>::diffusion(Time t, const Array_t<T> &x) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(x.size() == 1, "1-D array required");
#endif
    Matrix_t<T> m(1, 1, diffusion(t, x[0]));
    return m;
}

template <class T>
inline Disposable<Array_t<T> >
StochasticProcess1D_t<T>::expectation(Time t0, const Array_t<T> &x0,
                                      Time dt) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(x0.size() == 1, "1-D array required");
#endif
    Array_t<T> a(1, expectation(t0, x0[0], dt));
    return a;
}

template <class T>
inline Disposable<Matrix_t<T> >
StochasticProcess1D_t<T>::stdDeviation(Time t0, const Array_t<T> &x0,
                                       Time dt) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(x0.size() == 1, "1-D array required");
#endif
    Matrix_t<T> m(1, 1, stdDeviation(t0, x0[0], dt));
    return m;
}

template <class T>
inline Disposable<Matrix_t<T> >
StochasticProcess1D_t<T>::covariance(Time t0, const Array_t<T> &x0,
                                     Time dt) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(x0.size() == 1, "1-D array required");
#endif
    Matrix_t<T> m(1, 1, variance(t0, x0[0], dt));
    return m;
}

template <class T>
inline Disposable<Array_t<T> >
StochasticProcess1D_t<T>::evolve(Time t0, const Array_t<T> &x0, Time dt,
                                 const Array_t<T> &dw) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(x0.size() == 1, "1-D array required");
    QL_REQUIRE(dw.size() == 1, "1-D array required");
#endif
    Array_t<T> a(1, evolve(t0, x0[0], dt, dw[0]));
    return a;
}

template <class T>
inline Disposable<Array_t<T> >
StochasticProcess1D_t<T>::apply(const Array_t<T> &x0,
                                const Array_t<T> &dx) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(x0.size() == 1, "1-D array required");
    QL_REQUIRE(dx.size() == 1, "1-D array required");
#endif
    Array_t<T> a(1, apply(x0[0], dx[0]));
    return a;
}

// implementation

// base class

template <class T> StochasticProcess_t<T>::StochasticProcess_t() {}

template <class T>
StochasticProcess_t<T>::StochasticProcess_t(
    const boost::shared_ptr<discretization> &disc)
    : discretization_(disc) {}

template <class T> Size StochasticProcess_t<T>::factors() const {
    return size();
}

template <class T>
Disposable<Array_t<T> >
StochasticProcess_t<T>::expectation(Time t0, const Array_t<T> &x0,
                                    Time dt) const {
    return apply(x0, discretization_->drift(*this, t0, x0, dt));
}

template <class T>
Disposable<Matrix_t<T> >
StochasticProcess_t<T>::stdDeviation(Time t0, const Array_t<T> &x0,
                                     Time dt) const {
    return discretization_->diffusion(*this, t0, x0, dt);
}

template <class T>
Disposable<Matrix_t<T> >
StochasticProcess_t<T>::covariance(Time t0, const Array_t<T> &x0,
                                   Time dt) const {
    return discretization_->covariance(*this, t0, x0, dt);
}

template <class T>
Disposable<Array_t<T> >
StochasticProcess_t<T>::evolve(Time t0, const Array_t<T> &x0, Time dt,
                               const Array_t<T> &dw) const {
    return apply(expectation(t0, x0, dt), stdDeviation(t0, x0, dt) * dw);
}

template <class T>
Disposable<Array_t<T> >
StochasticProcess_t<T>::apply(const Array_t<T> &x0,
                              const Array_t<T> &dx) const {
    return x0 + dx;
}

template <class T> Time StochasticProcess_t<T>::time(const Date &) const {
    QL_FAIL("date/time conversion not supported");
}

template <class T> void StochasticProcess_t<T>::update() { notifyObservers(); }

// 1-D specialization

template <class T> StochasticProcess1D_t<T>::StochasticProcess1D_t() {}

template <class T>
StochasticProcess1D_t<T>::StochasticProcess1D_t(
    const boost::shared_ptr<discretization> &disc)
    : discretization_(disc) {}

template <class T>
T StochasticProcess1D_t<T>::expectation(Time t0, T x0, Time dt) const {
    return apply(x0, discretization_->drift(*this, t0, x0, dt));
}

template <class T>
T StochasticProcess1D_t<T>::stdDeviation(Time t0, T x0, Time dt) const {
    return discretization_->diffusion(*this, t0, x0, dt);
}

template <class T>
T StochasticProcess1D_t<T>::variance(Time t0, T x0, Time dt) const {
    return discretization_->variance(*this, t0, x0, dt);
}

template <class T>
T StochasticProcess1D_t<T>::evolve(Time t0, T x0, Time dt, T dw) const {
    return apply(expectation(t0, x0, dt), stdDeviation(t0, x0, dt) * dw);
}

template <class T> T StochasticProcess1D_t<T>::apply(Real x0, T dx) const {
    return x0 + dx;
}
}

#endif
