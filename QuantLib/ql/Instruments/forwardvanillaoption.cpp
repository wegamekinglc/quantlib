
/*
 Copyright (C) 2002, 2003 Ferdinando Ametrano

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

/*! \file forwardvanillaoption.cpp
    \brief Forward version of a vanilla option
*/

#include <ql/Instruments/forwardvanillaoption.hpp>

namespace QuantLib {

    namespace Instruments {

        ForwardVanillaOption::ForwardVanillaOption(Option::Type type,
            const RelinkableHandle<MarketElement>& underlying,
            const RelinkableHandle<TermStructure>& dividendTS,
            const RelinkableHandle<TermStructure>& riskFreeTS,
            const Exercise& exercise,
            const RelinkableHandle<BlackVolTermStructure>& volTS,
            const Handle<PricingEngine>& engine,
            double moneyness,
            Date resetDate,
            const std::string& isinCode,
            const std::string& description)
        : VanillaOption(type, underlying, 0.0, dividendTS, riskFreeTS,
          exercise, volTS, engine, isinCode, description),
          moneyness_(moneyness), resetDate_(resetDate) {}

        void ForwardVanillaOption::setupArguments(Arguments* args) const {
            VanillaOption::setupArguments(args);
            ForwardVanillaOption::arguments* arguments =
                dynamic_cast<ForwardVanillaOption::arguments*>(args);
            QL_REQUIRE(arguments != 0,
                       "ForwardVanillaOption::setupArguments :"
                       "wrong argument type");

            arguments->moneyness = moneyness_;
            arguments->resetDate = resetDate_;

        }

        void ForwardVanillaOption::performCalculations() const {
            if (isExpired()) {
                NPV_ = delta_ = gamma_ = theta_ =
                    vega_ =   rho_ = dividendRho_ = strikeSensitivity_ = 0.0;
            } else {
                Option::performCalculations();

                const ForwardVanillaOption::results* results =
                    dynamic_cast<const ForwardVanillaOption::results*>(
                        engine_->results());
                QL_ENSURE(results != 0,
                          "no results returned from pricing engine");
                delta_       = results->delta;
                gamma_       = results->gamma;
                theta_       = results->theta;
                vega_        = results->vega;
                rho_         = results->rho;
                dividendRho_ = results->dividendRho;

            }
            QL_ENSURE(NPV_ != Null<double>(),
                      "null value returned from option pricer");
        }

    }

}

