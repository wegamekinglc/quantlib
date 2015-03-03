/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

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

/*! \file bisection.hpp
    \brief bisection 1-D solver
*/

#ifndef quantlib_solver1d_bisection_h
#define quantlib_solver1d_bisection_h

#include <ql/math/solver1d.hpp>

namespace QuantLib {

    //! %Bisection 1-D solver
    /*! \test the correctness of the returned values is tested by
              checking them against known good results.
    */
    template<class T>
    class Bisection_t : public Solver1D<Bisection_t<T>,T> {
      public:
        template <class F>
        T solveImpl(const F& f,
                       T xAccuracy) const {

            /* The implementation of the algorithm was inspired by
               Press, Teukolsky, Vetterling, and Flannery,
               "Numerical Recipes in C", 2nd edition, Cambridge
               University Press
            */

            T dx, xMid, fMid;

            // Orient the search so that f>0 lies at root_+dx
            if (this->fxMin_ < 0.0) {
                dx = this->xMax_-this->xMin_;
                this->root_ = this->xMin_;
            } else {
                dx = this->xMin_-this->xMax_;
                this->root_ = this->xMax_;
            }

            while (this->evaluationNumber_<=this->maxEvaluations_) {
                dx /= 2.0;
                xMid = this->root_+dx;
                fMid = f(xMid);
                ++this->evaluationNumber_;
                if (fMid <= 0.0)
                    this->root_ = xMid;
                if (QLFCT::abs(dx) < xAccuracy || (close(fMid, T(0.0)))) {
                    f(this->root_);
                    ++this->evaluationNumber_;
                    return this->root_;
                }
            }
            QL_FAIL("maximum number of function evaluations ("
                    << this->maxEvaluations_ << ") exceeded");
        }
    };

    typedef Bisection_t<Real> Bisection;

}

#endif
