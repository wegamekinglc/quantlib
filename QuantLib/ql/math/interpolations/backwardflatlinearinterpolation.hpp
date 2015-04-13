/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

/*! \file bilinearinterpolation.hpp
    \brief backflat interpolation in first component, linear in second component
*/

#ifndef quantlib_backwardflatlinear_interpolation_hpp
#define quantlib_backwardflatlinear_interpolation_hpp

#include <ql/math/interpolations/interpolation2d.hpp>

namespace QuantLib {

namespace detail {

template <class I1, class I2, template <class> class M, class T>
class BackwardflatLinearInterpolationImpl_t
    : public Interpolation2D_t<T>::template templateImpl<I1, I2, M> {
  public:
    BackwardflatLinearInterpolationImpl_t(const I1 &xBegin, const I1 &xEnd,
                                          const I2 &yBegin, const I2 &yEnd,
                                          const M<T> &zData)
        : Interpolation2D_t<T>::template templateImpl<I1, I2, M>(xBegin, xEnd, yBegin, yEnd,
                                                   zData) {
        calculate();
    }
    void calculate() {}
    T value(T x, T y) const {
        Size j = this->locateY(y);
        T z1, z2;
        if (x <= this->xBegin_[0]) {
            z1 = this->zData_[j][0];
            z2 = this->zData_[j + 1][0];
        } else {
            Size i = this->locateX(x);
            if (x == this->xBegin_[i]) {
                z1 = this->zData_[j][i];
                z2 = this->zData_[j + 1][i];
            } else {
                z1 = this->zData_[j][i + 1];
                z2 = this->zData_[j + 1][i + 1];
            }
        }

        T u =
            (y - this->yBegin_[j]) / (this->yBegin_[j + 1] - this->yBegin_[j]);

        return (1.0 - u) * z1 + u * z2;
    }
};
}

template <class T = Real>
class BackwardflatLinearInterpolation_t : public Interpolation2D_t<T> {
  public:
    /*! \pre the \f$ x \f$ and \f$ y \f$ values must be sorted. */
    template <class I1, class I2, template <class> class M>
    BackwardflatLinearInterpolation_t(const I1 &xBegin, const I1 &xEnd,
                                      const I2 &yBegin, const I2 &yEnd,
                                      const M<T> &zData) {
        this->impl_ = boost::shared_ptr<typename Interpolation2D_t<T>::Impl>(
            new detail::BackwardflatLinearInterpolationImpl_t<I1, I2, M, T>(
                xBegin, xEnd, yBegin, yEnd, zData));
    }
};

typedef BackwardflatLinearInterpolation_t<Real> BackwardflatLinearInterpolation;

template <class T = Real> class BackwardflatLinear {
  public:
    template <class I1, class I2, template <class> class M>
    Interpolation2D_t<T> interpolate(const I1 &xBegin, const I1 &xEnd,
                                     const I2 &yBegin, const I2 &yEnd,
                                     const M<T> &z) const {
        return BackwardflatLinearInterpolation_t<T>(xBegin, xEnd, yBegin, yEnd, z);
    }
};
}

#endif
