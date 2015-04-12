/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003, 2006 Ferdinando Ametrano
 Copyright (C) 2004, 2005, 2006, 2007 StatPro Italia srl
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

/*! \file interpolation2d.hpp
    \brief abstract base classes for 2-D interpolations
*/

#ifndef quantlib_interpolation2D_hpp
#define quantlib_interpolation2D_hpp

#include <ql/math/interpolations/extrapolation.hpp>
#include <ql/math/comparison.hpp>
#include <ql/math/matrix.hpp>
#include <ql/errors.hpp>
#include <ql/types.hpp>
#include <vector>

namespace QuantLib {

//! base class for 2-D interpolations.
/*! Classes derived from this class will provide interpolated
    values from two sequences of length \f$ N \f$ and \f$ M \f$,
    representing the discretized values of the \f$ x \f$ and \f$ y
    \f$ variables, and a \f$ N \times M \f$ matrix representing
    the tabulated function values.
*/
template <class T> class Interpolation2D_t : public Extrapolator {
  protected:
    //! abstract base class for 2-D interpolation implementations
    class Impl {
      public:
        virtual ~Impl() {}
        virtual void calculate() = 0;
        virtual T xMin() const = 0;
        virtual T xMax() const = 0;
        virtual std::vector<T> xValues() const = 0;
        virtual Size locateX(T x) const = 0;
        virtual T yMin() const = 0;
        virtual T yMax() const = 0;
        virtual std::vector<T> yValues() const = 0;
        virtual Size locateY(T y) const = 0;
        virtual const Matrix_t<T> &zData() const = 0;
        virtual bool isInRange(T x, T y) const = 0;
        virtual T value(T x, T y) const = 0;
    };
    boost::shared_ptr<Impl> impl_;

  public:
    typedef T first_argument_type;
    typedef T second_argument_type;
    typedef T result_type;
    //! basic template implementation
    template <class I1, class I2, template <class> class M>
    class templateImpl : public Impl {
      public:
        templateImpl(const I1 &xBegin, const I1 &xEnd, const I2 &yBegin,
                     const I2 &yEnd, const M<T> &zData)
            : xBegin_(xBegin), xEnd_(xEnd), yBegin_(yBegin), yEnd_(yEnd),
              zData_(zData) {
            QL_REQUIRE(xEnd_ - xBegin_ >= 2,
                       "not enough x points to interpolate: at least 2 "
                       "required, "
                           << xEnd_ - xBegin_ << " provided");
            QL_REQUIRE(yEnd_ - yBegin_ >= 2,
                       "not enough y points to interpolate: at least 2 "
                       "required, "
                           << yEnd_ - yBegin_ << " provided");
        }
        T xMin() const { return *xBegin_; }
        T xMax() const { return *(xEnd_ - 1); }
        std::vector<T> xValues() const {
            return std::vector<T>(xBegin_, xEnd_);
        }
        T yMin() const { return *yBegin_; }
        T yMax() const { return *(yEnd_ - 1); }
        std::vector<T> yValues() const {
            return std::vector<T>(yBegin_, yEnd_);
        }
        const Matrix_t<T> &zData() const { return zData_; }
        bool isInRange(T x, T y) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
            for (I1 i = xBegin_, j = xBegin_ + 1; j != xEnd_; ++i, ++j)
                QL_REQUIRE(*j > *i, "unsorted x values");
#endif
            T x1 = xMin(), x2 = xMax();
            bool xIsInrange =
                (x >= x1 && x <= x2) || close(x, x1) || close(x, x2);
            if (!xIsInrange)
                return false;

#if defined(QL_EXTRA_SAFETY_CHECKS)
            for (I2 k = yBegin_, l = yBegin_ + 1; l != yEnd_; ++k, ++l)
                QL_REQUIRE(*l > *k, "unsorted y values");
#endif
            T y1 = yMin(), y2 = yMax();
            return (y >= y1 && y <= y2) || close(y, y1) || close(y, y2);
        }

      protected:
        Size locateX(T x) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
            for (I1 i = xBegin_, j = xBegin_ + 1; j != xEnd_; ++i, ++j)
                QL_REQUIRE(*j > *i, "unsorted x values");
#endif
            if (x < *xBegin_)
                return 0;
            else if (x > *(xEnd_ - 1))
                return xEnd_ - xBegin_ - 2;
            else
                return std::upper_bound(xBegin_, xEnd_ - 1, x) - xBegin_ - 1;
        }
        Size locateY(T y) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
            for (I2 k = yBegin_, l = yBegin_ + 1; l != yEnd_; ++k, ++l)
                QL_REQUIRE(*l > *k, "unsorted y values");
#endif
            if (y < *yBegin_)
                return 0;
            else if (y > *(yEnd_ - 1))
                return yEnd_ - yBegin_ - 2;
            else
                return std::upper_bound(yBegin_, yEnd_ - 1, y) - yBegin_ - 1;
        }
        I1 xBegin_, xEnd_;
        I2 yBegin_, yEnd_;
        const M<T> &zData_;
    };

  public:
    Interpolation2D_t() {}
    T operator()(T x, T y, bool allowExtrapolation = false) const {
        checkRange(x, y, allowExtrapolation);
        return impl_->value(x, y);
    }
    T xMin() const { return impl_->xMin(); }
    T xMax() const { return impl_->xMax(); }
    std::vector<T> xValues() const { return impl_->xValues(); }
    Size locateX(T x) const { return impl_->locateX(x); }
    T yMin() const { return impl_->yMin(); }
    T yMax() const { return impl_->yMax(); }
    std::vector<T> yValues() const { return impl_->yValues(); }
    Size locateY(T y) const { return impl_->locateY(y); }
    const Matrix_t<T> &zData() const { return impl_->zData(); }
    bool isInRange(T x, T y) const { return impl_->isInRange(x, y); }
    void update() { impl_->calculate(); }

  protected:
    void checkRange(T x, T y, bool extrapolate) const {
        QL_REQUIRE(extrapolate || allowsExtrapolation() ||
                       impl_->isInRange(x, y),
                   "interpolation range is ["
                       << impl_->xMin() << ", " << impl_->xMax() << "] x ["
                       << impl_->yMin() << ", " << impl_->yMax()
                       << "]: extrapolation at (" << x << ", " << y
                       << ") not allowed");
    }
};

typedef Interpolation2D_t<Real> Interpolation2D;
}

#endif
