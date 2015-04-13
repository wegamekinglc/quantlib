/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2006 Cristina Duminuco
 Copyright (C) 2006 Marco Bianchetti
 Copyright (C) 2007 StatPro Italia srl
 Copyright (C) 2014 Ferdinando Ametrano
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

/*! \file swaption.hpp
    \brief Swaption class
*/

#ifndef quantlib_instruments_swaption_base_hpp
#define quantlib_instruments_swaption_base_hpp

#include <ql/option.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>

namespace QuantLib {

class PricingEngine;

//! %settlement information
struct Settlement {
    enum Type { Physical, Cash };
};

std::ostream &operator<<(std::ostream &out, Settlement::Type type);

//! %Swaption class
/*! \ingroup instruments

    \test
    - the correctness of the returned value is tested by checking
      that the price of a payer (resp. receiver) swaption
      decreases (resp. increases) with the strike.
    - the correctness of the returned value is tested by checking
      that the price of a payer (resp. receiver) swaption
      increases (resp. decreases) with the spread.
    - the correctness of the returned value is tested by checking
      it against that of a swaption on a swap with no spread and a
      correspondingly adjusted fixed rate.
    - the correctness of the returned value is tested by checking
      it against a known good value.
    - the correctness of the returned value of cash settled swaptions
      is tested by checking the modified annuity against a value
      calculated without using the Swaption class.


    \todo add greeks and explicit exercise lag
*/
template <class T> class Swaption_t : public Option_t<T> {
  public:
    class arguments;
    class engine;
    Swaption_t(const boost::shared_ptr<VanillaSwap_t<T> > &swap,
               const boost::shared_ptr<Exercise> &exercise,
               Settlement::Type delivery = Settlement::Physical);
    //! \name Instrument interface
    //@{
    bool isExpired() const;
    void setupArguments(PricingEngine::arguments *) const;
    //@}
    //! \name Inspectors
    //@{
    Settlement::Type settlementType() const { return settlementType_; }
    typename VanillaSwap_t<T>::Type type() const { return swap_->type(); }
    const boost::shared_ptr<VanillaSwap_t<T> > &underlyingSwap() const {
        return swap_;
    }
    //@}
    //! implied volatility
    T impliedVolatility(T price,
                        const Handle<YieldTermStructure_t<T> > &discountCurve,
                        T guess, T accuracy = 1.0e-4,
                        Natural maxEvaluations = 100, T minVol = 1.0e-7,
                        T maxVol = 4.0, T displacement = 0.0) const;

  private:
    // arguments
    boost::shared_ptr<VanillaSwap_t<T> > swap_;
    // Handle<YieldTermStructure_t<T>> termStructure_;
    Settlement::Type settlementType_;
};

typedef Swaption_t<Real> Swaption;

//! %Arguments for swaption calculation
template <class T>
class Swaption_t<T>::arguments : public VanillaSwap_t<T>::arguments,
                                 public Option_t<T>::arguments {
  public:
    arguments() : settlementType(Settlement::Physical) {}
    boost::shared_ptr<VanillaSwap_t<T> > swap;
    Settlement::Type settlementType;
    void validate() const;
};

//! base class for swaption engines
template <class T>
class Swaption_t<T>::engine
    : public GenericEngine<typename Swaption_t<T>::arguments,
                           typename Swaption_t<T>::results> {};

} // namespace QuantLib

#endif
