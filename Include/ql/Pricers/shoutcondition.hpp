
/*
 * Copyright (C) 2000-2001 QuantLib Group
 *
 * This file is part of QuantLib.
 * QuantLib is a C++ open source library for financial quantitative
 * analysts and developers --- http://quantlib.org/
 *
 * QuantLib is free software and you are allowed to use, copy, modify, merge,
 * publish, distribute, and/or sell copies of it under the conditions stated
 * in the QuantLib License.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
 *
 * You should have received a copy of the license along with this file;
 * if not, contact ferdinando@ametrano.net
 * The license is also available at http://quantlib.org/LICENSE.TXT
 *
 * The members of the QuantLib Group are listed in the Authors.txt file, also
 * available at http://quantlib.org/Authors.txt
*/

/*! \file shoutcondition.hpp
    \brief shout option exercise condition

    $Id$
*/

// $Source$
// $Log$
// Revision 1.8  2001/07/25 15:47:28  sigmud
// Change from quantlib.sourceforge.net to quantlib.org
//
// Revision 1.7  2001/07/05 12:35:09  enri
// - added some static_cast<int>() to prevent gcc warnings
// - added some virtual constructor (same reason)
//
// Revision 1.6  2001/06/22 16:38:15  lballabio
// Improved documentation
//
// Revision 1.5  2001/05/25 09:29:40  nando
// smoothing #include xx.hpp and cutting old Log messages
//
// Revision 1.4  2001/05/24 15:38:08  nando
// smoothing #include xx.hpp and cutting old Log messages
//

#ifndef quantlib_pricers_shout_condition_h
#define quantlib_pricers_shout_condition_h

#include "ql/FiniteDifferences/standardstepcondition.hpp"

namespace QuantLib {

    namespace Pricers {

        class ShoutCondition 
        : public FiniteDifferences::StandardStepCondition {
          public:
            ShoutCondition(const Array& initialPrices, Time resTime,
                           Rate rate);
            void applyTo(Array& a, Time t) const;
          private:
            Rate rate_;
            Time resTime_;
            Array initialPrices_;
        };

        inline ShoutCondition::ShoutCondition(
            const Array& initialPrices, Time resTime, Rate rate)
            : rate_(rate), resTime_(resTime), initialPrices_(initialPrices){}

        inline void ShoutCondition::applyTo(Array& a, Time t) const {
            for (int i = 0; i < a.size(); i++)
                a[i] = QL_MAX(a[i], QL_EXP(-rate_ * (t - resTime_)) * 
                                           initialPrices_[i] );
        }

    }

}


#endif
