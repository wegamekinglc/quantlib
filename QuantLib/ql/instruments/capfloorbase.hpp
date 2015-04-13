/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2006, 2014 Ferdinando Ametrano
 Copyright (C) 2006 François du Vignaud
 Copyright (C) 2006, 2007 StatPro Italia srl
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

/*! \file capfloor.hpp
    \brief cap and floor class
*/

#ifndef quantlib_instruments_capfloor_base_hpp
#define quantlib_instruments_capfloor_base_hpp

#include <ql/instrument.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/handle.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/utilities/dataformatters.hpp>

namespace QuantLib {

//! Base class for cap-like instruments
/*! \ingroup instruments

    \test
    - the correctness of the returned value is tested by checking
      that the price of a cap (resp. floor) decreases
      (resp. increases) with the strike rate.
    - the relationship between the values of caps, floors and the
      resulting collars is checked.
    - the put-call parity between the values of caps, floors and
      swaps is checked.
    - the correctness of the returned implied volatility is tested
      by using it for reproducing the target value.
    - the correctness of the returned value is tested by checking
      it against a known good value.
*/
template <class T> class CapFloor_t : public Instrument {
  public:
    enum Type { Cap, Floor, Collar };
    class arguments;
    class engine;
    CapFloor_t(Type type, const typename Leg_t<T>::Type &floatingLeg,
               const std::vector<T> &capRates,
               const std::vector<T> &floorRates);
    CapFloor_t(Type type, const typename Leg_t<T>::Type &floatingLeg,
               const std::vector<T> &strikes);
    //! \name Instrument interface
    //@{
    bool isExpired() const;
    void setupArguments(PricingEngine::arguments *) const;
    //@}
    //! \name Inspectors
    //@{
    Type type() const { return type_; }
    const std::vector<T> &capRates() const { return capRates_; }
    const std::vector<T> &floorRates() const { return floorRates_; }
    const typename Leg_t<T>::Type &floatingLeg() const { return floatingLeg_; }

    Date startDate() const;
    Date maturityDate() const;
    boost::shared_ptr<FloatingRateCoupon_t<T> > lastFloatingRateCoupon() const;
    //! Returns the n-th optionlet as a new CapFloor with only one cash flow.
    boost::shared_ptr<CapFloor_t<T> > optionlet(const Size n) const;
    //@}
    T atmRate(const YieldTermStructure_t<T> &discountCurve) const;
    //! implied term volatility
    T impliedVolatility(T price, const Handle<YieldTermStructure_t<T> > &disc,
                        T guess, T accuracy = 1.0e-4,
                        Natural maxEvaluations = 100, T minVol = 1.0e-7,
                        T maxVol = 4.0, T displacement = 0.0) const;

  private:
    Type type_;
    typename Leg_t<T>::Type floatingLeg_;
    std::vector<T> capRates_;
    std::vector<T> floorRates_;
};

typedef CapFloor_t<Real> CapFloor;

namespace {
template <class T> class ImpliedVolHelperCapFloor {
  public:
    ImpliedVolHelperCapFloor(
        const CapFloor_t<T> &,
        const Handle<YieldTermStructure_t<T> > &discountCurve, T targetValue,
        T displacement);
    T operator()(T x) const;
    T derivative(T x) const;

  private:
    boost::shared_ptr<PricingEngine> engine_;
    Handle<YieldTermStructure_t<T> > discountCurve_;
    T targetValue_;
    boost::shared_ptr<SimpleQuote_t<T> > vol_;
    const Instrument::results *results_;
};
}

//! Concrete cap class
/*! \ingroup instruments */
template <class T> class Cap_t : public CapFloor_t<T> {
  public:
    Cap_t(const typename Leg_t<T>::Type &floatingLeg,
          const std::vector<T> &exerciseRates)
        : CapFloor_t<T>(CapFloor_t<T>::Cap, floatingLeg, exerciseRates,
                        std::vector<T>()) {}
};

typedef Cap_t<Real> Cap;

//! Concrete floor class
/*! \ingroup instruments */
template <class T> class Floor_t : public CapFloor_t<T> {
  public:
    Floor_t(const typename Leg_t<T>::Type &floatingLeg,
            const std::vector<T> &exerciseRates)
        : CapFloor_t<T>(CapFloor_t<T>::Floor, floatingLeg, std::vector<T>(),
                        exerciseRates) {}
};

typedef Floor_t<Real> Floor;

//! Concrete collar class
/*! \ingroup instruments */
template <class T> class Collar_t : public CapFloor_t<T> {
  public:
    Collar_t(const typename Leg_t<T>::Type &floatingLeg,
             const std::vector<T> &capRates, const std::vector<T> &floorRates)
        : CapFloor_t<T>(CapFloor_t<T>::Collar, floatingLeg, capRates,
                        floorRates) {}
};

typedef Collar_t<Real> Collar;

//! %Arguments for cap/floor calculation
template <class T>
class CapFloor_t<T>::arguments : public virtual PricingEngine::arguments {
  public:
    arguments() : type(CapFloor_t<T>::Type(-1)) {}
    CapFloor_t<T>::Type type;
    std::vector<Date> startDates;
    std::vector<Date> fixingDates;
    std::vector<Date> endDates;
    std::vector<Time> accrualTimes;
    std::vector<T> capRates;
    std::vector<T> floorRates;
    std::vector<T> forwards;
    std::vector<T> gearings;
    std::vector<T> spreads;
    std::vector<T> nominals;
    std::vector<boost::shared_ptr<InterestRateIndex_t<T> > > indexes;
    void validate() const;
};

//! base class for cap/floor engines
template <class T>
class CapFloor_t<T>::engine
    : public GenericEngine<CapFloor_t<T>::arguments, CapFloor_t<T>::results> {};

template <class T>
std::ostream &operator<<(std::ostream &, typename CapFloor_t<T>::Type);

// implementation

template <class T>
CapFloor_t<T>::CapFloor_t(CapFloor_t<T>::Type type,
                          const typename Leg_t<T>::Type &floatingLeg,
                          const std::vector<T> &capRates,
                          const std::vector<T> &floorRates)
    : type_(type), floatingLeg_(floatingLeg), capRates_(capRates),
      floorRates_(floorRates) {
    if (type_ == Cap || type_ == Collar) {
        QL_REQUIRE(!capRates_.empty(), "no cap rates given");
        capRates_.reserve(floatingLeg_.size());
        while (capRates_.size() < floatingLeg_.size())
            capRates_.push_back(capRates_.back());
    }
    if (type_ == Floor || type_ == Collar) {
        QL_REQUIRE(!floorRates_.empty(), "no floor rates given");
        floorRates_.reserve(floatingLeg_.size());
        while (floorRates_.size() < floatingLeg_.size())
            floorRates_.push_back(floorRates_.back());
    }
    typename Leg_t<T>::Type::const_iterator i;
    for (i = floatingLeg_.begin(); i != floatingLeg_.end(); ++i)
        this->registerWith(*i);

    this->registerWith(Settings::instance().evaluationDate());
}

template <class T>
CapFloor_t<T>::CapFloor_t(CapFloor_t<T>::Type type,
                          const typename Leg_t<T>::Type &floatingLeg,
                          const std::vector<T> &strikes)
    : type_(type), floatingLeg_(floatingLeg) {
    QL_REQUIRE(!strikes.empty(), "no strikes given");
    if (type_ == Cap) {
        capRates_ = strikes;
        capRates_.reserve(floatingLeg_.size());
        while (capRates_.size() < floatingLeg_.size())
            capRates_.push_back(capRates_.back());
    } else if (type_ == Floor) {
        floorRates_ = strikes;
        floorRates_.reserve(floatingLeg_.size());
        while (floorRates_.size() < floatingLeg_.size())
            floorRates_.push_back(floorRates_.back());
    } else
        QL_FAIL("only Cap/Floor types allowed in this constructor");

    typename Leg_t<T>::Type::const_iterator i;
    for (i = floatingLeg_.begin(); i != floatingLeg_.end(); ++i)
        this->registerWith(*i);

    this->registerWith(Settings::instance().evaluationDate());
}

template <class T> bool CapFloor_t<T>::isExpired() const {
    for (Size i = floatingLeg_.size(); i > 0; --i)
        if (!floatingLeg_[i - 1]->hasOccurred())
            return false;
    return true;
}

template <class T> Date CapFloor_t<T>::startDate() const {
    return CashFlows::startDate<T>(floatingLeg_);
}

template <class T> Date CapFloor_t<T>::maturityDate() const {
    return CashFlows::maturityDate<T>(floatingLeg_);
}

template <class T>
shared_ptr<FloatingRateCoupon_t<T> >
CapFloor_t<T>::lastFloatingRateCoupon() const {
    shared_ptr<CashFlow> lastCF(floatingLeg_.back());
    shared_ptr<FloatingRateCoupon_t<T> > lastFloatingCoupon =
        boost::dynamic_pointer_cast<FloatingRateCoupon_t<T> >(lastCF);
    return lastFloatingCoupon;
}

template <class T>
shared_ptr<CapFloor_t<T> > CapFloor_t<T>::optionlet(const Size i) const {
    QL_REQUIRE(i < floatingLeg().size(),
               io::ordinal(i + 1) << " optionlet does not exist, only "
                                  << floatingLeg().size());
    typename Leg_t<T>::Type cf(1, floatingLeg()[i]);

    std::vector<T> cap, floor;
    if (type() == Cap || type() == Collar)
        cap.push_back(capRates()[i]);
    if (type() == Floor || type() == Collar)
        floor.push_back(floorRates()[i]);

    return shared_ptr<CapFloor_t<T> >(
        new CapFloor_t<T>(type(), cf, cap, floor));
}

template <class T>
void CapFloor_t<T>::setupArguments(PricingEngine::arguments *args) const {
    CapFloor_t<T>::arguments *arguments =
        dynamic_cast<CapFloor_t<T>::arguments *>(args);
    QL_REQUIRE(arguments != 0, "wrong argument type");

    Size n = floatingLeg_.size();

    arguments->startDates.resize(n);
    arguments->fixingDates.resize(n);
    arguments->endDates.resize(n);
    arguments->accrualTimes.resize(n);
    arguments->forwards.resize(n);
    arguments->nominals.resize(n);
    arguments->gearings.resize(n);
    arguments->capRates.resize(n);
    arguments->floorRates.resize(n);
    arguments->spreads.resize(n);
    arguments->indexes.resize(n);

    arguments->type = type_;

    Date today = Settings::instance().evaluationDate();

    for (Size i = 0; i < n; ++i) {
        shared_ptr<FloatingRateCoupon_t<T> > coupon =
            boost::dynamic_pointer_cast<FloatingRateCoupon_t<T> >(
                floatingLeg_[i]);
        QL_REQUIRE(coupon, "non-FloatingRateCoupon given");
        arguments->startDates[i] = coupon->accrualStartDate();
        arguments->fixingDates[i] = coupon->fixingDate();
        arguments->endDates[i] = coupon->date();

        // this is passed explicitly for precision
        arguments->accrualTimes[i] = coupon->accrualPeriod();

        // this is passed explicitly for precision...
        if (arguments->endDates[i] >= today) { // ...but only if needed
            arguments->forwards[i] = coupon->adjustedFixing();
        } else {
            arguments->forwards[i] = Null<T>();
        }

        arguments->nominals[i] = coupon->nominal();
        T spread = coupon->spread();
        T gearing = coupon->gearing();
        arguments->gearings[i] = gearing;
        arguments->spreads[i] = spread;

        if (type_ == Cap || type_ == Collar)
            arguments->capRates[i] = (capRates_[i] - spread) / gearing;
        else
            arguments->capRates[i] = Null<T>();

        if (type_ == Floor || type_ == Collar)
            arguments->floorRates[i] = (floorRates_[i] - spread) / gearing;
        else
            arguments->floorRates[i] = Null<T>();

        arguments->indexes[i] = coupon->index();
    }
}

template <class T> void CapFloor_t<T>::arguments::validate() const {
    QL_REQUIRE(endDates.size() == startDates.size(),
               "number of start dates ("
                   << startDates.size()
                   << ") different from that of end dates (" << endDates.size()
                   << ")");
    QL_REQUIRE(accrualTimes.size() == startDates.size(),
               "number of start dates ("
                   << startDates.size()
                   << ") different from that of accrual times ("
                   << accrualTimes.size() << ")");
    QL_REQUIRE(
        type == CapFloor_t<T>::Floor || capRates.size() == startDates.size(),
        "number of start dates (" << startDates.size()
                                  << ") different from that of cap rates ("
                                  << capRates.size() << ")");
    QL_REQUIRE(
        type == CapFloor_t<T>::Cap || floorRates.size() == startDates.size(),
        "number of start dates (" << startDates.size()
                                  << ") different from that of floor rates ("
                                  << floorRates.size() << ")");
    QL_REQUIRE(gearings.size() == startDates.size(),
               "number of start dates ("
                   << startDates.size() << ") different from that of gearings ("
                   << gearings.size() << ")");
    QL_REQUIRE(spreads.size() == startDates.size(),
               "number of start dates (" << startDates.size()
                                         << ") different from that of spreads ("
                                         << spreads.size() << ")");
    QL_REQUIRE(nominals.size() == startDates.size(),
               "number of start dates ("
                   << startDates.size() << ") different from that of nominals ("
                   << nominals.size() << ")");
    QL_REQUIRE(forwards.size() == startDates.size(),
               "number of start dates ("
                   << startDates.size() << ") different from that of forwards ("
                   << forwards.size() << ")");
}

template <class T>
T CapFloor_t<T>::atmRate(const YieldTermStructure_t<T> &discountCurve) const {
    bool includeSettlementDateFlows = false;
    Date settlementDate = discountCurve.referenceDate();
    return CashFlows::atmRate(floatingLeg_, discountCurve,
                              includeSettlementDateFlows, settlementDate);
}

template <class T>
std::ostream &operator<<(std::ostream &out, typename CapFloor_t<T>::Type t) {
    switch (t) {
    case CapFloor_t<T>::Cap:
        return out << "Cap";
    case CapFloor_t<T>::Floor:
        return out << "Floor";
    case CapFloor_t<T>::Collar:
        return out << "Collar";
    default:
        QL_FAIL("unknown CapFloor_t<T>::Type (" << Integer(t) << ")");
    }
}

template <class T>
T CapFloor_t<T>::impliedVolatility(T targetValue,
                                   const Handle<YieldTermStructure_t<T> > &d,
                                   T guess, T accuracy, Natural maxEvaluations,
                                   T minVol, T maxVol, T displacement) const {
    // calculate();
    QL_REQUIRE(!isExpired(), "instrument expired");

    ImpliedVolHelperCapFloor<T> f(*this, d, targetValue, displacement);
    // Brent solver;
    NewtonSafe_t<T> solver;
    solver.setMaxEvaluations(maxEvaluations);
    return solver.solve(f, accuracy, guess, minVol, maxVol);
}

} // namespace QuantLib

#endif
