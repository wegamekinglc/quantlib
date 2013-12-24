/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013 Jose Aparicio

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

/*! \file blackcmcdsengine.hpp
    \brief Black constant maturity credit default swap engine
*/

#ifndef quantlib_black_cmcds_engine_hpp
#define quantlib_black_cmcds_engine_hpp

#include <ql/experimental/credit/constantmaturitycds.hpp>

namespace QuantLib {

    class Quote;
    class DefaultProbabilityTermStructure;

    /* Prices npv and default past events within the instrument protection
	 period, this embeds jump to default computation in the pricing.
    */
    class BalckConstantMaturityCDSEngine : public ConstantMaturityCDS::engine {
    public:
        BalckConstantMaturityCDSEngine(
              const Handle<DefaultProbabilityTermStructure>&,
              Real recoveryRate,
              const Handle<YieldTermStructure>& discountCurve,
              const Handle<Quote>& vol,
              boost::optional<bool> includeSettlementDateFlows = boost::none);
        void calculate() const;
      private:
        Handle<DefaultProbabilityTermStructure> probability_;
        Real recoveryRate_;
        Handle<Quote> vol_;
        Handle<YieldTermStructure> discountCurve_;
        boost::optional<bool> includeSettlementDateFlows_;
    };

}

#endif
