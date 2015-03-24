/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

/*! \file
\brief abstract base classes for 2-D flat extrapolations
*/

#ifndef quantlib_flatextrapolation2D_hpp
#define quantlib_flatextrapolation2D_hpp

#include <ql/math/interpolations/interpolation2d.hpp>

namespace QuantLib {

template <class T> class FlatExtrapolator2D_t : public Interpolation2D_t<T> {
  public:
    FlatExtrapolator2D_t(
        boost::shared_ptr<Interpolation2D_t<T> > decoratedInterpolation) {
        this->impl_ = boost::shared_ptr<typename Interpolation2D_t<T>::Impl>(
            new FlatExtrapolator2DImpl_t(decoratedInterpolation));
    }

  protected:
    class FlatExtrapolator2DImpl_t : public Interpolation2D_t<T>::Impl {
      public:
        FlatExtrapolator2DImpl_t(
            boost::shared_ptr<Interpolation2D_t<T> > decoratedInterpolation)
            : decoratedInterp_(decoratedInterpolation) {
            calculate();
        }
        T xMin() const { return decoratedInterp_->xMin(); }
        T xMax() const { return decoratedInterp_->xMax(); }
        std::vector<T> xValues() const { return decoratedInterp_->xValues(); }
        Size locateX(T x) const { return decoratedInterp_->locateX(x); }
        T yMin() const { return decoratedInterp_->yMin(); }
        T yMax() const { return decoratedInterp_->yMax(); }
        std::vector<T> yValues() const { return decoratedInterp_->yValues(); }
        Size locateY(T y) const { return decoratedInterp_->locateY(y); }
        const Matrix &zData() const { return decoratedInterp_->zData(); }
        bool isInRange(T x, T y) const {
            return decoratedInterp_->isInRange(x, y);
        }
        void update() { decoratedInterp_->update(); }
        void calculate() {}
        T value(T x, T y) const {
            x = bindX(x);
            y = bindY(y);
            return decoratedInterp_->operator()(x, y);
        }

      private:
        boost::shared_ptr<Interpolation2D_t<T> > decoratedInterp_;

        T bindX(T x) const {
            if (x < xMin())
                return xMin();
            if (x > xMax())
                return xMax();
            return x;
        }
        T bindY(T y) const {
            if (y < yMin())
                return yMin();
            if (y > yMax())
                return yMax();
            return y;
        }
    };
};

typedef FlatExtrapolator2D_t<Real> FlatExtrapolator2D;
}

#endif
