
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it under the
 terms of the QuantLib license.  You should have received a copy of the
 license along with this program; if not, please email quantlib-dev@lists.sf.net
 The license is also available online at http://quantlib.org/html/license.html

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file singleassetoption.cpp
    \brief common code for option evaluation
*/

#include <ql/Pricers/singleassetoption.hpp>
#include <ql/Solvers1D/brent.hpp>

namespace QuantLib {

    namespace Pricers {

        const double SingleAssetOption::dVolMultiplier_ = 0.0001;
        const double SingleAssetOption::dRMultiplier_   = 0.0001;
//        const double SingleAssetOption::dSMultiplier_   = 0.0001;
//        const double SingleAssetOption::dTMultiplier_   = 0.0001;

        SingleAssetOption::SingleAssetOption(Option::Type type,
            double underlying, double strike, Spread dividendYield,
            Rate riskFreeRate, Time residualTime, double volatility)
            : underlying_(underlying),
            payoff_(type, strike), dividendYield_(dividendYield),
            residualTime_(residualTime), hasBeenCalculated_(false),
            rhoComputed_(false), dividendRhoComputed_(false),
            vegaComputed_(false), thetaComputed_(false) {
            QL_REQUIRE(strike > 0.0,
                "SingleAssetOption::SingleAssetOption : strike ("+
                 DoubleFormatter::toString(strike)+
                 ") must be positive");
            QL_REQUIRE(underlying > 0.0,
                "SingleAssetOption::SingleAssetOption : underlying ("+
                 DoubleFormatter::toString(underlying)+
                 ") must be positive");
            QL_REQUIRE(residualTime > 0.0,
                "SingleAssetOption::SingleAssetOption : residual time ("+
                 DoubleFormatter::toString(residualTime)+
                 ") must be positive");
            //! Checks on volatility values are in setVolatility
            setVolatility(volatility);
            //! Checks on the risk-free rate are in setRiskFreeRate
            setRiskFreeRate(riskFreeRate);
        }

        void SingleAssetOption::setVolatility(double volatility) {
            QL_REQUIRE(volatility >= QL_MIN_VOLATILITY,
                 "SingleAssetOption::setVolatility ("+
                 DoubleFormatter::toString(volatility)+
                 "): Volatility too small");

            QL_REQUIRE(volatility <= QL_MAX_VOLATILITY,
                "SingleAssetOption::setVolatility ("+
                 DoubleFormatter::toString(volatility)+
                 ") : Volatility too high");

            volatility_ = volatility;
            hasBeenCalculated_ = false;
            vegaComputed_ = false;
            rhoComputed_ = false;
            dividendRhoComputed_ = false;
            thetaComputed_ = false;
        }

        void SingleAssetOption::setRiskFreeRate(Rate newRiskFreeRate) {
            riskFreeRate_ = newRiskFreeRate;
            hasBeenCalculated_ = false;
        }

        void SingleAssetOption::setDividendYield(Rate newDividendYield) {
            dividendYield_ = newDividendYield;
            hasBeenCalculated_ = false;
        }

        double SingleAssetOption::theta() const {

            if(!thetaComputed_) {

                // use Black-Scholes equation for theta computation
                theta_ =  riskFreeRate_ * value()
                        -(riskFreeRate_ - dividendYield_ ) * underlying_ * delta()
                        - 0.5 * volatility_ * volatility_ *
                                underlying_ * underlying_ * gamma();
                thetaComputed_ = true;
            }
            return theta_;
        }

        double SingleAssetOption::vega() const {

            if(!vegaComputed_) {

                double valuePlus = value();

                Handle<SingleAssetOption> brandNewFD = clone();
                double volMinus = volatility_ * (1.0 - dVolMultiplier_);
                brandNewFD -> setVolatility(volMinus);
                double valueMinus = brandNewFD -> value();

                vega_ = (valuePlus - valueMinus )/
                        (volatility_ * dVolMultiplier_);
                vegaComputed_ = true;
            }
            return vega_;
        }

        double SingleAssetOption::dividendRho() const {

            if(!dividendRhoComputed_){
                double valuePlus = value();

                Handle<SingleAssetOption> brandNewFD = clone();
                Rate dMinus = (dividendYield_ ?
                    dividendYield_ * (1.0 - dRMultiplier_) : 0.0001);
                brandNewFD -> setDividendYield(dMinus);
                double valueMinus = brandNewFD -> value();

                dividendRho_=(valuePlus - valueMinus) /
                    (dividendYield_ - dMinus);
                dividendRhoComputed_ = true;
            }
            return dividendRho_;
        }

        double SingleAssetOption::rho() const {

            if(!rhoComputed_){
                double valuePlus = value();

                Handle<SingleAssetOption> brandNewFD = clone();
                Rate rMinus= (riskFreeRate_ ?
                    riskFreeRate_ * (1.0 - dRMultiplier_) : 0.0001);
                brandNewFD -> setRiskFreeRate(rMinus);
                double valueMinus = brandNewFD -> value();

                rho_=(valuePlus - valueMinus) /
                    (riskFreeRate_ - rMinus);
                rhoComputed_  = true;
            }
            return rho_;
        }

        double SingleAssetOption::impliedVolatility(double targetValue, 
                    double accuracy, Size maxEvaluations, double minVol, 
                    double maxVol) const {
            // check option targetValue boundary condition
            QL_REQUIRE(targetValue > 0.0,
             "SingleAssetOption::impliedVol : targetValue must be positive");
            double optionValue = value();
            if (optionValue == targetValue)
                return volatility_;
            // clone used for root finding
            Handle<SingleAssetOption> tempBSM = clone();
            // objective function
            VolatilityFunction bsmf(tempBSM, targetValue);
            // solver
            Solvers1D::Brent s1d = Solvers1D::Brent();
            s1d.setMaxEvaluations(maxEvaluations);
            s1d.setLowerBound(minVol);
            s1d.setUpperBound(maxVol);

            return s1d.solve(bsmf, accuracy, volatility_, minVol, maxVol);
        }

        double SingleAssetOption::impliedDivYield(double targetValue, 
                    double accuracy, Size maxEvaluations, double minDivYield, 
                    double maxDivYield) const {
            // check option targetValue boundary condition
            QL_REQUIRE(targetValue > 0.0,
             "SingleAssetOption::impliedYield : targetValue must be positive");
            double optionValue = value();
            if (optionValue == targetValue)
                return dividendYield_;
            // clone used for root finding
            Handle<SingleAssetOption> tempBSM = clone();
            // objective function
            DivYieldFunction bsmf(tempBSM, targetValue);
            // solver
            Solvers1D::Brent s1d = Solvers1D::Brent();
            s1d.setMaxEvaluations(maxEvaluations);
            s1d.setLowerBound(minDivYield);
            s1d.setUpperBound(maxDivYield);    

            return s1d.solve(bsmf, accuracy, dividendYield_, 
                             minDivYield, maxDivYield);
        }

    }

}
