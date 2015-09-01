/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

#include <ql/experimental/fx/proxyfxtarfengine.hpp>
#include <ql/math/matrix.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/math/interpolations/bilinearinterpolation.hpp>
#include <ql/math/interpolations/flatextrapolation2d.hpp>

namespace QuantLib {

void ProxyFxTarfEngine::calculate() const {

    // handle the trivial cases
    FxTarfEngine::calculate();

    // determine the number of open fixings
    Date today = Settings::instance().evaluationDate();

    Size numberOpenFixings =
        std::distance(std::upper_bound(arguments_.openFixingDates.begin(),
                                       arguments_.openFixingDates.end(), today),
                      arguments_.openFixingDates.end());

    // determine the accumulated amount index
    Size accInd = std::upper_bound(proxy_->accBucketLimits.begin(),
                                   proxy_->accBucketLimits.end(),
                                   arguments_.accumulatedAmount) -
                  proxy_->accBucketLimits.begin() - 1;

    // sanity checks
    QL_REQUIRE(
        today >= proxy_->origEvalDate,
        "evaluation date ("
            << today
            << ") must be greater or equal than original evaluation date ("
            << proxy_->origEvalDate);
    QL_REQUIRE(numberOpenFixings <= proxy_->openFixingDates.size(),
               "number of open fixings ("
                   << numberOpenFixings
                   << ") must be less or equal the number of open fixings "
                      "provided by the proxy object ("
                   << proxy_->openFixingDates.size() << ")");
    QL_REQUIRE(accInd >= 0 && accInd < proxy_->accBucketLimits.size(),
               "accumulated amount index ("
                   << accInd << ") out of range given by the proxy (0..."
                   << proxy_->accBucketLimits.size() << ")");

    // logging
    // std::cerr << "proxy engine: use function (openFixingsIndex, accIndex) =
    // ("
    //           << (numberOpenFixings - 1) << "," << accInd
    //           << ") discount(last payment date=" << proxy_->lastPaymentDate
    //           << ") is " << discount_->discount(proxy_->lastPaymentDate)
    //           << std::endl;

    // set the core (trusted) region as addtional result
    // todo in case of interpolation we should return the intersection of
    // the core regions of the adjacent nodes
    results_.additionalResults["coreRegionMin"] =
        proxy_->functions[numberOpenFixings - 1][accInd]->coreRegion().first;
    results_.additionalResults["coreRegionMax"] =
        proxy_->functions[numberOpenFixings - 1][accInd]->coreRegion().second;

    // get the proxy function and return the npv, on forward basis
    if (!interpolate_)
        results_.value =
            proxy_->functions[numberOpenFixings - 1][accInd]->operator()(
                exchangeRate_->value()) *
                discount_->discount(proxy_->lastPaymentDate) +
            unsettledAmountNpv_;

    if (interpolate_) {
        std::vector<Real> accumulatedAmounts;
        for (Size i = 1; i <= proxy_->accBucketLimits.size(); ++i) {
            Real bucketMid = 0.5 * ((i < proxy_->accBucketLimits.size()
                                         ? proxy_->accBucketLimits[i]
                                         : arguments_.target) +
                                    proxy_->accBucketLimits[i - 1]);
            accumulatedAmounts.push_back(bucketMid);
        }
        std::vector<Real> proxies;
        for (Size j = 0; j < proxy_->accBucketLimits.size(); ++j) {
            proxies.push_back(
                proxy_->functions[numberOpenFixings - 1][j]->operator()(
                    exchangeRate_->value()));
        }
        boost::shared_ptr<LinearInterpolation> in =
            boost::make_shared<LinearInterpolation>(accumulatedAmounts.begin(),
                                                    accumulatedAmounts.end(),
                                                    proxies.begin());
        // for(Size i=0;i<data.columns();++i) {
        //     std::cerr << *(data.row_begin(0)+i) << "-";
        // }
        in->enableExtrapolation();
        results_.value = in->operator()(arguments_.accumulatedAmount) *
                             discount_->discount(proxy_->lastPaymentDate) +
                         unsettledAmountNpv_;
        // std::cerr << "interpolated amount (accumulatedAmount=" <<
        // arguments_.accumulatedAmount << " = " <<
        // in->operator()(arguments_.accumulatedAmount) << std::endl;
    }
} // calculate
} // namespace QuantLib
