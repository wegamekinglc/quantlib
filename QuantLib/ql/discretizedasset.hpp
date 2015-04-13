/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2004, 2005, 2006 StatPro Italia srl
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

/*! \file discretizedasset.hpp
    \brief Discretized asset classes
*/

#ifndef quantlib_discretized_asset_hpp
#define quantlib_discretized_asset_hpp

#include <ql/numericalmethod.hpp>
#include <ql/math/comparison.hpp>
#include <ql/exercise.hpp>

namespace QuantLib {

//! Discretized asset class used by numerical methods
template <class T> class DiscretizedAsset_t {
  public:
    DiscretizedAsset_t()
        : latestPreAdjustment_(QL_MAX_REAL),
          latestPostAdjustment_(QL_MAX_REAL) {}
    virtual ~DiscretizedAsset_t() {}

    //! \name inspectors
    //@{
    Time time() const { return time_; }
    Time &time() { return time_; }

    const Array_t<T> &values() const { return values_; }
    Array_t<T> &values() { return values_; }

    const boost::shared_ptr<Lattice_t<T>> &method() const { return method_; }
    //@}

    /*! \name High-level interface

        Users of discretized assets should use these methods in
        order to initialize, evolve and take the present value of
        the assets.  They call the corresponding methods in the
        Lattice interface, to which we refer for
        documentation.

        @{
    */
    void initialize(const boost::shared_ptr<Lattice_t<T>> &, Time t);
    void rollback(Time to);
    void partialRollback(Time to);
    T presentValue();
    //@}

    /*! \name Low-level interface

        These methods (that developers should override when
        deriving from DiscretizedAsset) are to be used by
        numerical methods and not directly by users, with the
        exception of adjustValues(), preAdjustValues() and
        postAdjustValues() that can be used together with
        partialRollback().

        @{
    */

    /*! This method should initialize the asset values to an Array
        of the given size and with values depending on the
        particular asset.
    */
    virtual void reset(Size size) = 0;

    /*! This method will be invoked after rollback and before any
        other asset (i.e., an option on this one) has any chance to
        look at the values. For instance, payments happening at times
        already spanned by the rollback will be added here.

        This method is not virtual; derived classes must override
        the protected preAdjustValuesImpl() method instead.
    */
    void preAdjustValues();

    /*! This method will be invoked after rollback and after any
        other asset had their chance to look at the values. For
        instance, payments happening at the present time (and therefore
        not included in an option to be exercised at this time) will be
        added here.

        This method is not virtual; derived classes must override
        the protected postAdjustValuesImpl() method instead.
    */
    void postAdjustValues();

    /*! This method performs both pre- and post-adjustment */
    void adjustValues() {
        preAdjustValues();
        postAdjustValues();
    }

    /*! This method returns the times at which the numerical
        method should stop while rolling back the asset. Typical
        examples include payment times, exercise times and such.

        \note The returned values are not guaranteed to be sorted.
    */
    virtual std::vector<Time> mandatoryTimes() const = 0;
    //@}
  protected:
    /*! This method checks whether the asset was rolled at the
        given time. */
    bool isOnTime(Time t) const;
    /*! This method performs the actual pre-adjustment */
    virtual void preAdjustValuesImpl() {}
    /*! This method performs the actual post-adjustment */
    virtual void postAdjustValuesImpl() {}

    Time time_;
    Time latestPreAdjustment_, latestPostAdjustment_;
    Array_t<T> values_;

  private:
    boost::shared_ptr<Lattice_t<T>> method_;
};

typedef DiscretizedAsset_t<Real> DiscretizedAsset;

//! Useful discretized discount bond asset
template <class T>
class DiscretizedDiscountBond_t : public DiscretizedAsset_t<T> {
  public:
    DiscretizedDiscountBond_t() {}
    void reset(Size size) { this->values_ = Array_t<T>(size, 1.0); }
    std::vector<Time> mandatoryTimes() const { return std::vector<Time>(); }
};

typedef DiscretizedDiscountBond_t<Real> DiscretizedDiscountBond;

//! Discretized option on a given asset
/*! \warning it is advised that derived classes take care of
             creating and initializing themselves an instance of
             the underlying.
*/
template <class T> class DiscretizedOption_t : public DiscretizedAsset_t<T> {
  public:
    DiscretizedOption_t(
        const boost::shared_ptr<DiscretizedAsset_t<T> > &underlying,
        Exercise::Type exerciseType, const std::vector<Time> &exerciseTimes)
        : underlying_(underlying), exerciseType_(exerciseType),
          exerciseTimes_(exerciseTimes) {}
    void reset(Size size);
    std::vector<Time> mandatoryTimes() const;

  protected:
    void postAdjustValuesImpl();
    void applyExerciseCondition();
    boost::shared_ptr<DiscretizedAsset_t<T> > underlying_;
    Exercise::Type exerciseType_;
    std::vector<Time> exerciseTimes_;
};

typedef DiscretizedOption_t<Real> DiscretizedOption;

// inline definitions

template <class T>
inline void DiscretizedAsset_t<T>::initialize(
    const boost::shared_ptr<Lattice_t<T> > &method, Time t) {
    method_ = method;
    method_->initialize(*this, t);
}

template <class T> inline void DiscretizedAsset_t<T>::rollback(Time to) {
    method_->rollback(*this, to);
}

template <class T> inline void DiscretizedAsset_t<T>::partialRollback(Time to) {
    method_->partialRollback(*this, to);
}

template <class T> inline T DiscretizedAsset_t<T>::presentValue() {
    return method_->presentValue(*this);
}

template <class T> inline void DiscretizedAsset_t<T>::preAdjustValues() {
    if (!close_enough(time(), latestPreAdjustment_)) {
        preAdjustValuesImpl();
        latestPreAdjustment_ = time();
    }
}

template <class T> inline void DiscretizedAsset_t<T>::postAdjustValues() {
    if (!close_enough(time(), latestPostAdjustment_)) {
        postAdjustValuesImpl();
        latestPostAdjustment_ = time();
    }
}

template <class T> inline bool DiscretizedAsset_t<T>::isOnTime(Time t) const {
    const TimeGrid &grid = method()->timeGrid();
    return close_enough(grid[grid.index(t)], time());
}

template <class T> inline void DiscretizedOption_t<T>::reset(Size size) {
    QL_REQUIRE(this->method() == underlying_->method(),
               "option and underlying were initialized on "
               "different methods");
    this->values_ = Array_t<T>(size, 0.0);
    this->adjustValues();
}

template <class T>
inline std::vector<Time> DiscretizedOption_t<T>::mandatoryTimes() const {
    std::vector<Time> times = underlying_->mandatoryTimes();
    // discard negative times...
    std::vector<Time>::const_iterator i =
        std::find_if(exerciseTimes_.begin(), exerciseTimes_.end(),
                     std::bind2nd(std::greater_equal<Time>(), 0.0));
    // and add the positive ones
    times.insert(times.end(), i, exerciseTimes_.end());
    return times;
}

template <class T>
inline void DiscretizedOption_t<T>::applyExerciseCondition() {
    for (Size i = 0; i < this->values_.size(); i++)
        this->values_[i] =
            QLFCT::max(underlying_->values()[i], this->values_[i]);
}

// implementation

template <class T> void DiscretizedOption_t<T>::postAdjustValuesImpl() {
    /* In the real world, with time flowing forward, first
       any payment is settled and only after options can be
       exercised. Here, with time flowing backward, options
       must be exercised before performing the adjustment.
    */
    underlying_->partialRollback(this->time());
    underlying_->preAdjustValues();
    Size i;
    switch (exerciseType_) {
    case Exercise::American:
        if (this->time_ >= exerciseTimes_[0] &&
            this->time_ <= exerciseTimes_[1])
            applyExerciseCondition();
        break;
    case Exercise::Bermudan:
    case Exercise::European:
        for (i = 0; i < exerciseTimes_.size(); i++) {
            Time t = exerciseTimes_[i];
            if (t >= 0.0 && this->isOnTime(t))
                applyExerciseCondition();
        }
        break;
    default:
        QL_FAIL("invalid exercise type");
    }
    underlying_->postAdjustValues();
}
} // namespace QuantLib

#endif
