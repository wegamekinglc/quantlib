
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

/*! \file geometricasianoption.hpp
    \brief geometric Asian option

    $Id$
*/

// $Source$
// $Log$
// Revision 1.7  2001/07/25 15:47:28  sigmud
// Change from quantlib.sourceforge.net to quantlib.org
//
// Revision 1.6  2001/07/19 16:40:11  lballabio
// Improved docs a bit
//
// Revision 1.5  2001/06/22 16:38:15  lballabio
// Improved documentation
//
// Revision 1.4  2001/05/24 15:38:08  nando
// smoothing #include xx.hpp and cutting old Log messages
//

#ifndef    quantlib_geometric_asian_option_pricer_h
#define    quantlib_geometric_asian_option_pricer_h

#include "ql/Pricers/bsmeuropeanoption.hpp"

namespace QuantLib {

    namespace Pricers {

        //! geometric Asian option
        class GeometricAsianOption : public BSMEuropeanOption    {
           public:
            GeometricAsianOption(Type type, double underlying, double    strike,
                Rate dividendYield, Rate riskFreeRate, Time    residualTime,
                double volatility);
            double vega() const;
            double rho() const;
            Handle<BSMOption> clone() const;
        };


        // inline definitions
        
        inline GeometricAsianOption::GeometricAsianOption(Type type,
            double underlying, double strike, Rate dividendYield,
            Rate riskFreeRate, Time residualTime, double volatility):
            BSMEuropeanOption(type, underlying, strike, dividendYield/2,
            riskFreeRate/2-volatility*volatility/12, residualTime,
            volatility/QL_SQRT(3)){}

        inline double GeometricAsianOption::rho() const{
            return BSMEuropeanOption::rho()/2;
        }

        inline double GeometricAsianOption::vega() const{
            return BSMEuropeanOption::vega()/QL_SQRT(3)
                -BSMEuropeanOption::rho()*volatility_*volatility_/4;
        }

        inline Handle<BSMOption> GeometricAsianOption::clone() const{
            return Handle<BSMOption>(new GeometricAsianOption(*this));
        }

    }

}


#endif
