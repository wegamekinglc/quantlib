
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

/*! \file himalayapathpricer.hpp
    \brief multipath pricer for European-type Himalaya option

    $Id$
*/

// $Source$
// $Log$
// Revision 1.9  2001/07/25 15:47:28  sigmud
// Change from quantlib.sourceforge.net to quantlib.org
//
// Revision 1.8  2001/07/19 16:40:10  lballabio
// Improved docs a bit
//
// Revision 1.7  2001/07/13 14:29:08  sigmud
// removed a few gcc compile warnings
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

#ifndef quantlib_himalaya_path_pricer_h
#define quantlib_himalaya_path_pricer_h

#include <vector>

#include "ql/MonteCarlo/multipathpricer.hpp"

namespace QuantLib {

    namespace MonteCarlo {

        //! multipath pricer for European-type Himalaya option
        /*! The payoff of an himalaya option is computed in the following way:
            given a basket of N assets, and N time periods, at end of
            each period the option who performed the best is added to the
            average and then discarded from the basket. At the end of the
            N periods the option pays the max between the strike and the
            average of the best performers.
        */
        class HimalayaPathPricer : public MultiPathPricer {
          public:
            HimalayaPathPricer():MultiPathPricer(){}
            HimalayaPathPricer(const Array &underlying,
                double strike, double discount);
            double value(const MultiPath &path) const;
          protected:
            Array underlying_;
            double strike_;
            double discount_;
        };

    }

}

#endif
