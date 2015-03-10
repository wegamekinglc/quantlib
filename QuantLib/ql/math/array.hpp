/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2009 StatPro Italia srl
 Copyright (C) 2004 Ferdinando Ametrano

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

/*! \file array.hpp
    \brief 1-D array used in linear algebra.
*/

#ifndef quantlib_array_hpp
#define quantlib_array_hpp

#include <ql/types.hpp>
#include <ql/errors.hpp>
#include <ql/utilities/disposable.hpp>
#include <ql/utilities/null.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/scoped_array.hpp>
#include <boost/type_traits.hpp>
#include <functional>
#include <numeric>
#include <vector>
#include <iomanip>

namespace QuantLib {

//! 1-D array used in linear algebra.
/*! This class implements the concept of vector as used in linear
    algebra.
    As such, it is <b>not</b> meant to be used as a container -
    <tt>std::vector</tt> should be used instead.

    \test construction of arrays is checked in a number of cases
*/
template <class T> class Array_t {
  public:
    //! \name Constructors, destructor, and assignment
    //@{
    //! creates the array with the given dimension
    explicit Array_t(Size size = 0);
    //! creates the array and fills it with <tt>value</tt>
    Array_t(Size size, T value);
    /*! \brief creates the array and fills it according to
        \f$ a_{0} = value, a_{i}=a_{i-1}+increment \f$
    */
    Array_t(Size size, T value, T increment);
    Array_t(const Array_t<T> &);
    Array_t(const Disposable<Array_t<T> > &);
    //! creates the array from an iterable sequence
    template <class ForwardIterator>
    Array_t(ForwardIterator begin, ForwardIterator end);

    Array_t<T> &operator=(const Array_t<T> &);
    Array_t<T> &operator=(const Disposable<Array_t<T> > &);
    bool operator==(const Array_t<T> &) const;
    bool operator!=(const Array_t<T> &) const;
    //@}
    /*! \name Vector algebra

        <tt>v += x</tt> and similar operation involving a scalar value
        are shortcuts for \f$ \forall i : v_i = v_i + x \f$

        <tt>v *= w</tt> and similar operation involving two vectors are
        shortcuts for \f$ \forall i : v_i = v_i \times w_i \f$

        \pre all arrays involved in an algebraic expression must have
        the same size.
    */
    //@{
    const Array_t<T> &operator+=(const Array_t<T> &);
    const Array_t<T> &operator+=(T);
    const Array_t<T> &operator-=(const Array_t<T> &);
    const Array_t<T> &operator-=(T);
    const Array_t<T> &operator*=(const Array_t<T> &);
    const Array_t<T> &operator*=(T);
    const Array_t<T> &operator/=(const Array_t<T> &);
    const Array_t<T> &operator/=(T);
    //@}
    //! \name Element access
    //@{
    //! read-only
    T operator[](Size) const;
    T at(Size) const;
    T front() const;
    T back() const;
    //! read-write
    T &operator[](Size);
    T &at(Size);
    T &front();
    T &back();
    //@}
    //! \name Inspectors
    //@{
    //! dimension of the array
    Size size() const;
    //! whether the array is empty
    bool empty() const;
    //@}
    typedef Size size_type;
    typedef T value_type;
    typedef T *iterator;
    typedef const T *const_iterator;
    typedef boost::reverse_iterator<iterator> reverse_iterator;
    typedef boost::reverse_iterator<const_iterator> const_reverse_iterator;
    //! \name Iterator access
    //@{
    const_iterator begin() const;
    iterator begin();
    const_iterator end() const;
    iterator end();
    const_reverse_iterator rbegin() const;
    reverse_iterator rbegin();
    const_reverse_iterator rend() const;
    reverse_iterator rend();
    //@}
    //! \name Utilities
    //@{
    void swap(Array_t<T> &); // never throws
                             //@}

  private:
    boost::scoped_array<T> data_;
    Size n_;
};

//! specialization of null template for this class
template <class T> class Null<Array_t<T> > {
  public:
    Null() {}
    operator Array_t<T>() const { return Array_t<T>(); }
};

/*! \relates Array */
template <class T> T DotProduct(const Array_t<T> &, const Array_t<T> &);

// unary operators
/*! \relates Array */
template <class T> const Disposable<Array_t<T> > operator+(const Array_t<T> &v);
/*! \relates Array */
template <class T> const Disposable<Array_t<T> > operator-(const Array_t<T> &v);

// binary operators
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator+(const Array_t<T> &, const Array_t<T> &);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator+(const Array_t<T> &, T);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator+(T, const Array_t<T> &);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator-(const Array_t<T> &, const Array_t<T> &);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator-(const Array_t<T> &, T);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator-(T, const Array_t<T> &);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator*(const Array_t<T> &, const Array_t<T> &);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator*(const Array_t<T> &, T);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator*(T, const Array_t<T> &);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator/(const Array_t<T> &, const Array_t<T> &);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator/(const Array_t<T> &, T);
/*! \relates Array */
template <class T>
const Disposable<Array_t<T> > operator/(T, const Array_t<T> &);

// math functions
/*! \relates Array */
template <class T> const Disposable<Array_t<T> > Abs(const Array_t<T> &);
/*! \relates Array */
template <class T> const Disposable<Array_t<T> > Sqrt(const Array_t<T> &);
/*! \relates Array */
template <class T> const Disposable<Array_t<T> > Log(const Array_t<T> &);
/*! \relates Array */
template <class T> const Disposable<Array_t<T> > Exp(const Array_t<T> &);
/*! \relates Array */
template <class T> const Disposable<Array_t<T> > Pow(const Array_t<T> &, T);

// utilities
/*! \relates Array */
template <class T> void swap(Array_t<T> &, Array_t<T> &);

// format
/*! \relates Array */
template <class T> std::ostream &operator<<(std::ostream &, const Array_t<T> &);

// inline definitions

template <class T>
inline Array_t<T>::Array_t(Size size)
    : data_(size ? new T[size] : (T *)(0)), n_(size) {}

template <class T>
inline Array_t<T>::Array_t(Size size, T value)
    : data_(size ? new T[size] : (T *)(0)), n_(size) {
    std::fill(begin(), end(), value);
}

template <class T>
inline Array_t<T>::Array_t(Size size, T value, T increment)
    : data_(size ? new T[size] : (T *)(0)), n_(size) {
    for (iterator i = begin(); i != end(); i++, value += increment)
        *i = value;
}

template <class T>
inline Array_t<T>::Array_t(const Array_t<T> &from)
    : data_(from.n_ ? new T[from.n_] : (T *)(0)), n_(from.n_) {
#if defined(QL_PATCH_MSVC) && defined(QL_DEBUG)
    if (n_)
#endif
        std::copy(from.begin(), from.end(), begin());
}

template <class T>
inline Array_t<T>::Array_t(const Disposable<Array_t<T> > &from)
    : data_((T *)(0)), n_(0) {
    swap(const_cast<Disposable<Array_t<T> > &>(from));
}

namespace detail {

template <class T, class I>
inline void _fill_array_(Array_t<T> &a, boost::scoped_array<T> &data_, Size &n_,
                         I begin, I end, const boost::true_type &) {
    // we got redirected here from a call like Array(3, 4)
    // because it matched the constructor below exactly with
    // ForwardIterator = int.  What we wanted was fill an
    // Array with a given value, which we do here.
    Size n = begin;
    T value = end;
    data_.reset(n ? new T[n] : (T *)(0));
    n_ = n;
    std::fill(a.begin(), a.end(), value);
}

template <class T, class I>
inline void _fill_array_(Array_t<T> &a, boost::scoped_array<T> &data_, Size &n_,
                         I begin, I end, const boost::false_type &) {
    // true iterators
    Size n = std::distance(begin, end);
    data_.reset(n ? new T[n] : (T *)(0));
    n_ = n;
#if defined(QL_PATCH_MSVC) && defined(QL_DEBUG)
    if (n_)
#endif
        std::copy(begin, end, a.begin());
}
}

template <class T>
template <class ForwardIterator>
inline Array_t<T>::Array_t(ForwardIterator begin, ForwardIterator end) {
    // Unfortunately, calls such as Array(3, 4) match this constructor.
    // We have to detect integral types and dispatch.
    detail::_fill_array_(*this, data_, n_, begin, end,
                         boost::is_integral<ForwardIterator>());
}

template <class T>
inline Array_t<T> &Array_t<T>::operator=(const Array_t<T> &from) {
    // strong guarantee
    Array_t<T> temp(from);
    swap(temp);
    return *this;
}

template <class T>
inline bool Array_t<T>::operator==(const Array_t<T> &to) const {
    return (n_ == to.n_) && std::equal(begin(), end(), to.begin());
}

template <class T>
inline bool Array_t<T>::operator!=(const Array_t<T> &to) const {
    return !(this->operator==(to));
}

template <class T>
inline Array_t<T> &Array_t<T>::operator=(const Disposable<Array_t<T> > &from) {
    swap(const_cast<Disposable<Array_t<T> > &>(from));
    return *this;
}

template <class T>
inline const Array_t<T> &Array_t<T>::operator+=(const Array_t<T> &v) {
    QL_REQUIRE(n_ == v.n_, "arrays with different sizes ("
                               << n_ << ", " << v.n_ << ") cannot be added");
    std::transform(begin(), end(), v.begin(), begin(), std::plus<T>());
    return *this;
}

template <class T> inline const Array_t<T> &Array_t<T>::operator+=(T x) {
    std::transform(begin(), end(), begin(), std::bind2nd(std::plus<T>(), x));
    return *this;
}

template <class T>
inline const Array_t<T> &Array_t<T>::operator-=(const Array_t<T> &v) {
    QL_REQUIRE(n_ == v.n_, "arrays with different sizes ("
                               << n_ << ", " << v.n_
                               << ") cannot be subtracted");
    std::transform(begin(), end(), v.begin(), begin(), std::minus<T>());
    return *this;
}

template <class T> inline const Array_t<T> &Array_t<T>::operator-=(T x) {
    std::transform(begin(), end(), begin(), std::bind2nd(std::minus<T>(), x));
    return *this;
}

template <class T>
inline const Array_t<T> &Array_t<T>::operator*=(const Array_t<T> &v) {
    QL_REQUIRE(n_ == v.n_, "arrays with different sizes ("
                               << n_ << ", " << v.n_
                               << ") cannot be multiplied");
    std::transform(begin(), end(), v.begin(), begin(), std::multiplies<T>());
    return *this;
}

template <class T> inline const Array_t<T> &Array_t<T>::operator*=(T x) {
    std::transform(begin(), end(), begin(),
                   std::bind2nd(std::multiplies<T>(), x));
    return *this;
}

template <class T>
inline const Array_t<T> &Array_t<T>::operator/=(const Array_t<T> &v) {
    QL_REQUIRE(n_ == v.n_, "arrays with different sizes ("
                               << n_ << ", " << v.n_ << ") cannot be divided");
    std::transform(begin(), end(), v.begin(), begin(), std::divides<T>());
    return *this;
}

template <class T> inline const Array_t<T> &Array_t<T>::operator/=(T x) {
    std::transform(begin(), end(), begin(), std::bind2nd(std::divides<T>(), x));
    return *this;
}

template <class T> inline T Array_t<T>::operator[](Size i) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < n_, "index (" << i << ") must be less than " << n_
                                 << ": array access out of range");
#endif
    return data_.get()[i];
}

template <class T> inline T Array_t<T>::at(Size i) const {
    QL_REQUIRE(i < n_, "index (" << i << ") must be less than " << n_
                                 << ": array access out of range");
    return data_.get()[i];
}

template <class T> inline T Array_t<T>::front() const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(n_ > 0, "null Array_t<T>: array access out of range");
#endif
    return data_.get()[0];
}

template <class T> inline T Array_t<T>::back() const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(n_ > 0, "null Array_t<T>: array access out of range");
#endif
    return data_.get()[n_ - 1];
}

template <class T> inline T &Array_t<T>::operator[](Size i) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < n_, "index (" << i << ") must be less than " << n_
                                 << ": array access out of range");
#endif
    return data_.get()[i];
}

template <class T> inline T &Array_t<T>::at(Size i) {
    QL_REQUIRE(i < n_, "index (" << i << ") must be less than " << n_
                                 << ": array access out of range");
    return data_.get()[i];
}

template <class T> inline T &Array_t<T>::front() {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(n_ > 0, "null Array: array access out of range");
#endif
    return data_.get()[0];
}

template <class T> inline T &Array_t<T>::back() {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(n_ > 0, "null Array: array access out of range");
#endif
    return data_.get()[n_ - 1];
}

template <class T> inline Size Array_t<T>::size() const { return n_; }

template <class T> inline bool Array_t<T>::empty() const { return n_ == 0; }

template <class T>
inline typename Array_t<T>::const_iterator Array_t<T>::begin() const {
    return data_.get();
}

template <class T> inline typename Array_t<T>::iterator Array_t<T>::begin() {
    return data_.get();
}

template <class T>
inline typename Array_t<T>::const_iterator Array_t<T>::end() const {
    return data_.get() + n_;
}

template <class T> inline typename Array_t<T>::iterator Array_t<T>::end() {
    return data_.get() + n_;
}

template <class T>
inline typename Array_t<T>::const_reverse_iterator Array_t<T>::rbegin() const {
    return const_reverse_iterator(end());
}

template <class T>
inline typename Array_t<T>::reverse_iterator Array_t<T>::rbegin() {
    return reverse_iterator(end());
}

template <class T>
inline typename Array_t<T>::const_reverse_iterator Array_t<T>::rend() const {
    return const_reverse_iterator(begin());
}

template <class T>
inline typename Array_t<T>::reverse_iterator Array_t<T>::rend() {
    return reverse_iterator(begin());
}

template <class T> inline void Array_t<T>::swap(Array_t<T> &from) {
    using std::swap;
    data_.swap(from.data_);
    swap(n_, from.n_);
}

// dot product

template <class T>
inline T DotProduct(const Array_t<T> &v1, const Array_t<T> &v2) {
    QL_REQUIRE(v1.size() == v2.size(), "arrays with different sizes ("
                                           << v1.size() << ", " << v2.size()
                                           << ") cannot be multiplied");
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

// overloaded operators

// unary

template <class T>
inline const Disposable<Array_t<T> > operator+(const Array_t<T> &v) {
    Array_t<T> result = v;
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator-(const Array_t<T> &v) {
    Array_t<T> result(v.size());
    std::transform(v.begin(), v.end(), result.begin(), std::negate<T>());
    return result;
}

// binary operators

template <class T>
inline const Disposable<Array_t<T> > operator+(const Array_t<T> &v1,
                                               const Array_t<T> &v2) {
    QL_REQUIRE(v1.size() == v2.size(), "arrays with different sizes ("
                                           << v1.size() << ", " << v2.size()
                                           << ") cannot be added");
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                   std::plus<T>());
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator+(const Array_t<T> &v1, T a) {
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), result.begin(),
                   std::bind2nd(std::plus<T>(), a));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator+(T a, const Array_t<T> &v2) {
    Array_t<T> result(v2.size());
    std::transform(v2.begin(), v2.end(), result.begin(),
                   std::bind1st(std::plus<T>(), a));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator-(const Array_t<T> &v1,
                                               const Array_t<T> &v2) {
    QL_REQUIRE(v1.size() == v2.size(), "arrays with different sizes ("
                                           << v1.size() << ", " << v2.size()
                                           << ") cannot be subtracted");
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                   std::minus<T>());
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator-(const Array_t<T> &v1, T a) {
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), result.begin(),
                   std::bind2nd(std::minus<T>(), a));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator-(T a, const Array_t<T> &v2) {
    Array_t<T> result(v2.size());
    std::transform(v2.begin(), v2.end(), result.begin(),
                   std::bind1st(std::minus<T>(), a));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator*(const Array_t<T> &v1,
                                               const Array_t<T> &v2) {
    QL_REQUIRE(v1.size() == v2.size(), "arrays with different sizes ("
                                           << v1.size() << ", " << v2.size()
                                           << ") cannot be multiplied");
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                   std::multiplies<T>());
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator*(const Array_t<T> &v1, T a) {
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), result.begin(),
                   std::bind2nd(std::multiplies<T>(), a));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator*(T a, const Array_t<T> &v2) {
    Array_t<T> result(v2.size());
    std::transform(v2.begin(), v2.end(), result.begin(),
                   std::bind1st(std::multiplies<T>(), a));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator/(const Array_t<T> &v1,
                                               const Array_t<T> &v2) {
    QL_REQUIRE(v1.size() == v2.size(), "arrays with different sizes ("
                                           << v1.size() << ", " << v2.size()
                                           << ") cannot be divided");
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                   std::divides<T>());
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator/(const Array_t<T> &v1, T a) {
    Array_t<T> result(v1.size());
    std::transform(v1.begin(), v1.end(), result.begin(),
                   std::bind2nd(std::divides<T>(), a));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator/(T a, const Array_t<T> &v2) {
    Array_t<T> result(v2.size());
    std::transform(v2.begin(), v2.end(), result.begin(),
                   std::bind1st(std::divides<T>(), a));
    return result;
}

// functions

template <class T>
inline const Disposable<Array_t<T> > Abs(const Array_t<T> &v) {
    Array_t<T> result(v.size());
    for(Size i=0;i<v.size();++i)
        result[i] = QLFCT::abs(v[i]);
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > Sqrt(const Array_t<T> &v) {
    Array_t<T> result(v.size());
    for(Size i=0;i<v.size();++i)
        result[i] = QLFCT::sqrt(v[i]);
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > Log(const Array_t<T> &v) {
    Array_t<T> result(v.size());
    for(Size i=0;i<v.size();++i)
        result[i] = QLFCT::log(v[i]);
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > Exp(const Array_t<T> &v) {
    Array_t<T> result(v.size());
    for(Size i=0;i<v.size();++i)
        result[i] = QLFCT::exp(v[i]);
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > Pow(const Array_t<T> &v, T alpha) {
    Array_t<T> result(v.size());
    for(Size i=0;i<v.size();++i)
        result[i] = QLFCT::pow(v[i],alpha);
    return result;
}

template <class T> inline void swap(Array_t<T> &v, Array_t<T> &w) { v.swap(w); }

template <class T>
inline std::ostream &operator<<(std::ostream &out, const Array_t<T> &a) {
    std::streamsize width = out.width();
    out << "[ ";
    if (!a.empty()) {
        for (Size n = 0; n < a.size() - 1; ++n)
            out << std::setw(int(width)) << a[n] << "; ";
        out << std::setw(int(width)) << a.back();
    }
    out << " ]";
    return out;
}

typedef Array_t<Real> Array;

} // namespace QuantLib

#endif
