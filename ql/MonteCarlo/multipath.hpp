
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

// $Id$

#ifndef quantlib_montecarlo_multi_path_h
#define quantlib_montecarlo_multi_path_h

#include <ql/MonteCarlo/path.hpp>
#include <vector>

namespace QuantLib {

    namespace MonteCarlo {

        //! MultiPath implements multiple factors evolving at the same time
        /*!

        MultiPath contains the list of variations for each asset,
        \f[
            \log \frac{Y^j_{i+1}}{Y^j_i} \mathrm{for} i = 0, \ldots, n-1
            \qquad \mathrm{and} \qquad j = 0, \ldots, m-1
        \f]
        where \f$ Y^j_i \f$ is the value of the underlying \f$ j \f$
        at discretized time \f$ t_i \f$. The first index refers to the
        underlying, the second to the time position MultiPath[j,i]

        \todo make it time-aware

        */

        //! single random walk
        class MultiPath {
          public:
            MultiPath(Size nAsset,
                      Size pathSize);
            MultiPath(const std::vector<Path>& multiPath);
            //! \name inspectors
            //@{
            Size assetNumber() const {return multiPath_.size(); }
            Size pathSize() const {return multiPath_[0].size(); }
            //@}
            //! \name read/write access to components
            //@{
            const Path& operator[](Size j) const {return multiPath_[j]; }
            Path& operator[](Size j) {return multiPath_[j]; }
            //@}
          private:
            std::vector<Path> multiPath_;
        };

        // inline definitions

        inline MultiPath::MultiPath(Size nAsset, Size pathSize)
            : multiPath_(nAsset,Path(pathSize)) {
            QL_REQUIRE(nAsset > 0,
                "MultiPath: number of asset must be > zero");
            QL_REQUIRE(pathSize > 0,
                "MultiPath: pathSize must be > zero");
        }

        inline MultiPath::MultiPath(const std::vector<Path>& multiPath)
            : multiPath_(multiPath) {}


    }

}


#endif
