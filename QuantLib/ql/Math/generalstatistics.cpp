
/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2003 RiskMap srl

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

/*! \file generalstatistics.cpp
    \brief statistics tool
*/

#include <ql/Math/generalstatistics.hpp>
#include <ql/Math/functional.hpp>

namespace QuantLib {

    namespace Math {

        double GeneralStatistics::weightSum() const {
            double result = 0.0;
            std::vector<std::pair<double,double> >::const_iterator it;
            for (it=samples_.begin(); it!=samples_.end(); it++) {
                result += it->second;
            }
            return result;
        }

        double GeneralStatistics::mean() const {
            Size N = samples();
            QL_ENSURE(N != 0,
                      "GeneralStatistics::mean() : "
                      "empty sample set");
            // eat our own dog food
            return expectationValue(identity<double>(),
                                    everywhere()).first;
        }

        double GeneralStatistics::variance() const {
            Size N = samples();
            QL_ENSURE(N > 1,
                      "GeneralStatistics::variance() : "
                      "sample number <=1, unsufficient");
            // Subtract the mean and square. Repeat on the whole range.
            // Hopefully, the whole thing will be inlined in a single loop.
            double s2 = expectationValue(compose(square<double>(),
                                                 std::bind2nd(
                                                     std::minus<double>(),
                                                     mean())),
                                         everywhere()).first;
            return s2*N/(N-1.0);
        }

        double GeneralStatistics::skewness() const {
            Size N = samples();
            QL_REQUIRE(N > 2,
                       "GeneralStatistics::skewness() : "
                       "sample number <=2, unsufficient");

            double x = expectationValue(compose(cube<double>(),
                                                std::bind2nd(
                                                    std::minus<double>(),
                                                    mean())),
                                        everywhere()).first;
            double sigma = standardDeviation();

            return (x/(sigma*sigma*sigma))*(N/(N-1.0))*(N/(N-2.0));
        }

        double GeneralStatistics::kurtosis() const {
            Size N = samples();
            QL_REQUIRE(N > 3,
                       "GeneralStatistics::kurtosis() : "
                       "sample number <=3, unsufficient");

            double x = expectationValue(compose(fourth_power<double>(),
                                                std::bind2nd(
                                                    std::minus<double>(),
                                                    mean())),
                                        everywhere()).first;
            double sigma2 = variance();

            double c1 = (N/(N-1.0)) * (N/(N-2.0)) * ((N+1.0)/(N-3.0));
            double c2 = 3.0 * ((N-1.0)/(N-2.0)) * ((N-1.0)/(N-3.0));

            return c1*(x/(sigma2*sigma2))-c2;
        }

        /*! \pre percent must be in range (0%-100%] */
        double GeneralStatistics::percentile(double percent) const {

            QL_REQUIRE(percent > 0.0 && percent <= 1.0,
                       "GeneralStatistics::percentile() : "
                       "percentile (" +
                        DoubleFormatter::toString(percent) +
                       ") must be in (0.0, 1.0]");

            double sampleWeight = weightSum();
            QL_REQUIRE(sampleWeight>0.0,
                       "GeneralStatistics::percentile() : "
                       "empty sample set");

            std::sort(samples_.begin(), samples_.end());

            std::vector<std::pair<double,double> >::iterator k, l;
            k = samples_.begin();
            l = samples_.end()-1;
            /* the sum of weight is non null, therefore there's
               at least one sample */
            double integral = k->second, target = percent*sampleWeight;
            while (integral < target && k != l) {
                k++;
                integral += k->second;
            }
            return k->first;

            /*
            // interpolating ... wrong thing to do, but here's the code
            if (k==samples_.begin()) {
                return k->first;
            } else {
                // just in case there are more samples at value k->first
                double lastAddedWeight = k->second;
                std::vector<std::pair<double,double> >::iterator kk = k;
                kk++;
                while (kk!=samples_.end() && kk->first==k->first) {
                    lastAddedWeight += kk->second;
                    integral        += kk->second;
                    kk++;
                }
                double lambda = (integral - perc) / lastAddedWeight;
                return (1.0-lambda) * (k->first) + lambda * ((k-1)->first);
            }
            */
        }

        /*! \pre percent must be in range (0%-100%] */
        double GeneralStatistics::topPercentile(double percent) const {

            QL_REQUIRE(percent > 0.0 && percent <= 1.0,
                       "GeneralStatistics::topPercentile() : "
                       "percentile (" +
                        DoubleFormatter::toString(percent) +
                       ") must be in (0.0, 1.0]");

            double sampleWeight = weightSum();
            QL_REQUIRE(sampleWeight > 0.0,
                       "GeneralStatistics::topPercentile() : "
                       "empty sample set");

            std::sort(samples_.begin(), samples_.end());

            std::vector<std::pair<double,double> >::reverse_iterator k, l;
            k = samples_.rbegin();
            l = samples_.rend()-1;
            /* the sum of weight is non null, therefore there's
               at least one sample */
            double integral = k->second, target = percent*sampleWeight;
            while (integral < target && k != l) {
                k++;
                integral += k->second;
            }
            return k->first;
        }

    }

}
