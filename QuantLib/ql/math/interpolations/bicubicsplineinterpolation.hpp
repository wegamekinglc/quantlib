/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2004 StatPro Italia srl
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

/*! \file bicubicsplineinterpolation.hpp
    \brief bicubic spline interpolation between discrete points
*/

#ifndef quantlib_bicubic_spline_interpolation_hpp
#define quantlib_bicubic_spline_interpolation_hpp

#include <ql/math/interpolations/interpolation2d.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>

namespace QuantLib {

namespace detail {

template <class T> class BicubicSplineDerivatives_t {
  public:
    virtual ~BicubicSplineDerivatives_t() {}
    virtual T derivativeX(T x, T y) const = 0;
    virtual T derivativeY(T x, T y) const = 0;
    virtual T derivativeXY(T x, T y) const = 0;
    virtual T secondDerivativeX(T x, T y) const = 0;
    virtual T secondDerivativeY(T x, T y) const = 0;
};

template <class I1, class I2, template <class> class M, class T = Real>
class BicubicSplineImpl_t : public Interpolation2D::templateImpl<I1, I2, M>,
                            public BicubicSplineDerivatives_t<T> {
  public:
    BicubicSplineImpl_t(const I1 &xBegin, const I1 &xEnd, const I2 &yBegin,
                        const I2 &yEnd, const M<T> &zData)
        : Interpolation2D_t<T>::template templateImpl<I1, I2, M>(
              xBegin, xEnd, yBegin, yEnd, zData) {
        calculate();
    }
    void calculate() {
        splines_.resize(this->zData_.rows());
        for (Size i = 0; i < (this->zData_.rows()); ++i)
            splines_[i] = CubicInterpolation_t<T>(
                this->xBegin_, this->xEnd_, this->zData_.row_begin(i),
                CubicInterpolation_t<T>::Spline, false,
                CubicInterpolation_t<T>::SecondDerivative, 0.0,
                CubicInterpolation_t<T>::SecondDerivative, 0.0);
    }
    T value(T x, T y) const {
        std::vector<T> section(splines_.size());
        for (Size i = 0; i < splines_.size(); i++)
            section[i] = splines_[i](x, true);

        CubicInterpolation_t<T> spline(
            this->yBegin_, this->yEnd_, section.begin(),
            CubicInterpolation_t<T>::Spline, false,
            CubicInterpolation_t<T>::SecondDerivative, 0.0,
            CubicInterpolation_t<T>::SecondDerivative, 0.0);
        return spline(y, true);
    }

    T derivativeX(T x, T y) const {
        std::vector<T> section(this->zData_.columns());
        for (Size i = 0; i < section.size(); ++i) {
            section[i] = value(this->xBegin_[i], y);
        }

        return CubicInterpolation_t<T>(
                   this->xBegin_, this->xEnd_, section.begin(),
                   CubicInterpolation_t<T>::Spline, false,
                   CubicInterpolation_t<T>::SecondDerivative, 0.0,
                   CubicInterpolation_t<T>::SecondDerivative,
                   0.0).derivative(x);
    }

    T secondDerivativeX(T x, T y) const {
        std::vector<T> section(this->zData_.columns());
        for (Size i = 0; i < section.size(); ++i) {
            section[i] = value(this->xBegin_[i], y);
        }

        return CubicInterpolation_t<T>(
                   this->xBegin_, this->xEnd_, section.begin(),
                   CubicInterpolation_t<T>::Spline, false,
                   CubicInterpolation_t<T>::SecondDerivative, 0.0,
                   CubicInterpolation_t<T>::SecondDerivative,
                   0.0).secondDerivative(x);
    }

    T derivativeY(T x, T y) const {
        std::vector<T> section(splines_.size());
        for (Size i = 0; i < splines_.size(); i++)
            section[i] = splines_[i](x, true);

        return CubicInterpolation_t<T>(
                   this->yBegin_, this->yEnd_, section.begin(),
                   CubicInterpolation_t<T>::Spline, false,
                   CubicInterpolation_t<T>::SecondDerivative, 0.0,
                   CubicInterpolation_t<T>::SecondDerivative,
                   0.0).derivative(y);
    }

    T secondDerivativeY(T x, T y) const {
        std::vector<T> section(splines_.size());
        for (Size i = 0; i < splines_.size(); i++)
            section[i] = splines_[i](x, true);

        return CubicInterpolation_t<T>(
                   this->yBegin_, this->yEnd_, section.begin(),
                   CubicInterpolation_t<T>::Spline, false,
                   CubicInterpolation_t<T>::SecondDerivative, 0.0,
                   CubicInterpolation_t<T>::SecondDerivative,
                   0.0).secondDerivative(y);
    }

    T derivativeXY(T x, T y) const {
        std::vector<T> section(this->zData_.columns());
        for (Size i = 0; i < section.size(); ++i) {
            section[i] = derivativeY(this->xBegin_[i], y);
        }

        return CubicInterpolation_t<T>(
                   this->xBegin_, this->xEnd_, section.begin(),
                   CubicInterpolation_t<T>::Spline, false,
                   CubicInterpolation_t<T>::SecondDerivative, 0.0,
                   CubicInterpolation_t<T>::SecondDerivative,
                   0.0).derivative(x);
    }

  private:
    std::vector<Interpolation_t<T> > splines_;
};
}

//! bicubic-spline interpolation between discrete points
/*! \todo revise end conditions */
template <class T> class BicubicSpline_t : public Interpolation2D_t<T> {
  public:
    /*! \pre the \f$ x \f$ and \f$ y \f$ values must be sorted. */
    template <class I1, class I2, template<class> class M>
    BicubicSpline_t(const I1 &xBegin, const I1 &xEnd, const I2 &yBegin,
                    const I2 &yEnd, const M<T> &zData) {
        this->impl_ = boost::shared_ptr<typename Interpolation2D_t<T>::Impl>(
            new detail::BicubicSplineImpl_t<I1, I2, M, T>(xBegin, xEnd, yBegin,
                                                          yEnd, zData));
    }

    T derivativeX(T x, T y) const {
        return boost::dynamic_pointer_cast<
                   detail::BicubicSplineDerivatives_t<T> >(this->impl_)
            ->derivativeX(x, y);
    }
    T derivativeY(T x, T y) const {
        return boost::dynamic_pointer_cast<
                   detail::BicubicSplineDerivatives_t<T> >(this->impl_)
            ->derivativeY(x, y);
    }
    T secondDerivativeX(T x, T y) const {
        return boost::dynamic_pointer_cast<
                   detail::BicubicSplineDerivatives_t<T> >(this->impl_)
            ->secondDerivativeX(x, y);
    }
    T secondDerivativeY(T x, T y) const {
        return boost::dynamic_pointer_cast<
                   detail::BicubicSplineDerivatives_t<T> >(this->impl_)
            ->secondDerivativeY(x, y);
    }

    T derivativeXY(T x, T y) const {
        return boost::dynamic_pointer_cast<
                   detail::BicubicSplineDerivatives_t<T> >(this->impl_)
            ->derivativeXY(x, y);
    }
};

typedef BicubicSpline_t<Real> BicubicSpline;

//! bicubic-spline-interpolation factory
template <class T = Real> class Bicubic {
  public:
    template <class I1, class I2, class M>
    Interpolation2D_t<T> interpolate(const I1 &xBegin, const I1 &xEnd,
                                     const I2 &yBegin, const I2 &yEnd,
                                     const M &z) const {
        return BicubicSpline_t<T>(xBegin, xEnd, yBegin, yEnd, z);
    }
};
}

#endif
