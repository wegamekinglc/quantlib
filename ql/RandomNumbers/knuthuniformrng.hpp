
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
/*! \file knuthuniformrng.hpp
    \brief Knuth uniform random number generator

    \fullpath
    ql/RandomNumbers/%knuthuniformrng.hpp
*/

// $Id$

#ifndef quantlib_knuth_uniform_rng_h
#define quantlib_knuth_uniform_rng_h

#include <ql/MonteCarlo/sample.hpp>
#include <vector>

namespace QuantLib {

    //! Random Number Generators and Low Discrepancy Sequences
    /*! See sect. \ref randomnumbers */
    namespace RandomNumbers {

        //! Uniform random number generator
        /*! Random number generator by Knuth.
            For more details see Knuth, Seminumerical Algorithms,
            3rd edition, Section 3.6.
            \note This is <b>not</b> Knuth's original implementation which
            is available at
            http://www-cs-faculty.stanford.edu/~knuth/programs.html,
            but rather a slightly modified version wrapped in a C++ class.
            Such modifications did not affect the code but only the data
            structures used, which were converted in their C++/STL
            equivalents.
        */
        class KnuthUniformRng {
          public:
            typedef MonteCarlo::Sample<double> sample_type;
            /*! if the given seed is 0, a random seed will be chosen
                based on clock() */
            explicit KnuthUniformRng(long seed = 0);
            /*! returns a sample with weight 1.0 containing a random number
                uniformly chosen from (0.0,1.0) */
            sample_type next() const;
          private:
            /* Knuth's names and routines were preserved as much as possible
               while changing the data structures to more modern ones. */
            static const int KK, LL, TT, QUALITY;
            mutable std::vector<double> ranf_arr_buf;
            mutable std::vector<double>::const_iterator ranf_arr_ptr,
                                                        ranf_arr_sentinel;
            mutable std::vector<double> ran_u;
            double mod_sum(double x, double y) const;
            bool is_odd(int s) const;
            void ranf_start(long seed);
            void ranf_array(std::vector<double>& aa, int n) const;
            double ranf_arr_cycle() const;
        };


        // inline definitions

        inline KnuthUniformRng::sample_type KnuthUniformRng::next() const {
            double result = (ranf_arr_ptr != ranf_arr_sentinel ?
                             *ranf_arr_ptr++ :
                             ranf_arr_cycle());
            return sample_type(result,1.0);
        }

        inline double KnuthUniformRng::mod_sum(double x, double y) const {
            return (x+y)-int(x+y);
        }

        inline bool KnuthUniformRng::is_odd(int s) const {
            return (s&1) != 0;
        }

    }

}


#endif
