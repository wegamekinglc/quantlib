
/*
 Copyright (C) 2000, 2001, 2002 RiskMap srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software developed by the QuantLib Group; you can
 redistribute it and/or modify it under the terms of the QuantLib License;
 either version 1.0, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 QuantLib License for more details.

 You should have received a copy of the QuantLib License along with this
 program; if not, please email ferdinando@ametrano.net

 The QuantLib License is also available at http://quantlib.org/license.html
 The members of the QuantLib Group are listed in the QuantLib License
*/
/*! \file floatingratecoupon.hpp
    \brief Coupon at par on a term structure

    \fullpath
    ql/CashFlows/%floatingratecoupon.hpp
*/

// $Id$

#ifndef quantlib_floating_rate_coupon_hpp
#define quantlib_floating_rate_coupon_hpp

#include <ql/CashFlows/coupon.hpp>
#include <ql/Indexes/xibor.hpp>

namespace QuantLib {

    namespace CashFlows {

        //! %coupon at par on a term structure
        /*! \warning This class does not perform any date adjustment,
            i.e., the start and end date passed upon construction
            should be already rolled to a business day.
        */
        class FloatingRateCoupon : public Coupon,
                                   public Patterns::Observer {
          public:
            FloatingRateCoupon(double nominal,
                const Handle<Indexes::Xibor>& index,
                const RelinkableHandle<TermStructure>& termStructure,
                const Date& startDate, const Date& endDate,
                int fixingDays,
                Spread spread = 0.0,
                const Date& refPeriodStart = Date(),
                const Date& refPeriodEnd = Date());
            ~FloatingRateCoupon();
            //! \name CashFlow interface
            //@{
            double amount() const;
            //@}
            //! \name Inspectors
            //@{
            const Handle<Indexes::Xibor>& index() const;
            Spread spread() const;
            //@}
            //! \name Observer interface
            //@{
            void update();
            //@}
          private:
            RelinkableHandle<TermStructure> termStructure_;
            Handle<Indexes::Xibor> index_;
            int fixingDays_;
            Spread spread_;
        };


        // inline definitions

        inline FloatingRateCoupon::~FloatingRateCoupon() {
            termStructure_.unregisterObserver(this);
        }

        inline const Handle<Indexes::Xibor>&
        FloatingRateCoupon::index() const {
            return index_;
        }

        inline Spread FloatingRateCoupon::spread() const {
            return spread_;
        }

        inline void FloatingRateCoupon::update() {
            notifyObservers();
        }

    }

}


#endif
