/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl
 Copyright (C) 2003, 2004 Ferdinando Ametrano
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

/*! \file matrix.hpp
    \brief matrix used in linear algebra.
*/

#ifndef quantlib_matrix_hpp
#define quantlib_matrix_hpp

#include <ql/math/array.hpp>
#include <ql/utilities/steppingiterator.hpp>

#if defined(QL_PATCH_MSVC)
#pragma warning(push)
#pragma warning(disable : 4180)
#pragma warning(disable : 4127)
#endif

#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

#if !defined(QL_NO_UBLAS_SUPPORT)
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#endif

#if defined(QL_PATCH_MSVC)
#pragma warning(pop)
#endif

#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic pop
#endif

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace QuantLib {

//! %Matrix used in linear algebra.
/*! This class implements the concept of Matrix as used in linear
    algebra. As such, it is <b>not</b> meant to be used as a
    container.
*/
template <class T> class Matrix_t {
  public:
    //! \name Constructors, destructor, and assignment
    //@{
    //! creates a null matrix
    Matrix_t();
    //! creates a matrix with the given dimensions
    Matrix_t(Size rows, Size columns);
    //! creates the matrix and fills it with <tt>value</tt>
    Matrix_t(Size rows, Size columns, T value);
    //! creates the matrix and fills with values given by iterator start and end
    template<class I> Matrix_t(Size rows, Size columns, I start, I end);
    Matrix_t(const Matrix_t<T> &);
    Matrix_t(const Disposable<Matrix_t<T> > &);
    Matrix_t<T> &operator=(const Matrix_t<T> &);
    Matrix_t<T> &operator=(const Disposable<Matrix_t<T> > &);
    //@}

    //! \name Algebraic operators
    /*! \pre all matrices involved in an algebraic expression must have
             the same size.
    */
    //@{
    const Matrix_t<T> &operator+=(const Matrix_t<T> &);
    const Matrix_t<T> &operator-=(const Matrix_t<T> &);
    const Matrix_t<T> &operator*=(T);
    const Matrix_t<T> &operator/=(T);
    //@}

    typedef T *iterator;
    typedef const T *const_iterator;
    typedef boost::reverse_iterator<iterator> reverse_iterator;
    typedef boost::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef T *row_iterator;
    typedef const T *const_row_iterator;
    typedef boost::reverse_iterator<row_iterator> reverse_row_iterator;
    typedef boost::reverse_iterator<const_row_iterator>
        const_reverse_row_iterator;
    typedef step_iterator<iterator> column_iterator;
    typedef step_iterator<const_iterator> const_column_iterator;
    typedef boost::reverse_iterator<column_iterator> reverse_column_iterator;
    typedef boost::reverse_iterator<const_column_iterator>
        const_reverse_column_iterator;
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
    const_row_iterator row_begin(Size i) const;
    row_iterator row_begin(Size i);
    const_row_iterator row_end(Size i) const;
    row_iterator row_end(Size i);
    const_reverse_row_iterator row_rbegin(Size i) const;
    reverse_row_iterator row_rbegin(Size i);
    const_reverse_row_iterator row_rend(Size i) const;
    reverse_row_iterator row_rend(Size i);
    const_column_iterator column_begin(Size i) const;
    column_iterator column_begin(Size i);
    const_column_iterator column_end(Size i) const;
    column_iterator column_end(Size i);
    const_reverse_column_iterator column_rbegin(Size i) const;
    reverse_column_iterator column_rbegin(Size i);
    const_reverse_column_iterator column_rend(Size i) const;
    reverse_column_iterator column_rend(Size i);
    //@}

    //! \name Element access
    //@{
    const_row_iterator operator[](Size) const;
    const_row_iterator at(Size) const;
    row_iterator operator[](Size);
    row_iterator at(Size);
    Disposable<Array_t<T> > diagonal(void) const;
    //@}

    //! \name Inspectors
    //@{
    Size rows() const;
    Size columns() const;
    bool empty() const;
    //@}

    //! \name Utilities
    //@{
    void swap(Matrix_t<T> &);
    //@}
  private:
    boost::scoped_array<T> data_;
    Size rows_, columns_;
};

// algebraic operators

/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > operator+(const Matrix_t<T> &,
                                         const Matrix_t<T> &);
/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > operator-(const Matrix_t<T> &,
                                         const Matrix_t<T> &);
/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > operator*(const Matrix_t<T> &, T);
/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > operator*(T, const Matrix_t<T> &);
/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > operator/(const Matrix_t<T> &, T);

// vectorial products

/*! \relates Matrix */
template <class T>
const Disposable<Array_t<T> > operator*(const Array_t<T> &,
                                        const Matrix_t<T> &);
/*! \relates Matrix */
template <class T>
const Disposable<Array_t<T> > operator*(const Matrix_t<T> &,
                                        const Array_t<T> &);
/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > operator*(const Matrix_t<T> &,
                                         const Matrix_t<T> &);

// misc. operations

/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > transpose(const Matrix_t<T> &);

/*! \relates Matrix */
template <class T>
const Disposable<Matrix_t<T> > outerProduct(const Array_t<T> &v1,
                                            const Array_t<T> &v2);

/*! \relates Matrix */
template <class T, class Iterator1, class Iterator2>
const Disposable<Matrix_t<T> > outerProduct(Iterator1 v1begin, Iterator1 v1end,
                                            Iterator2 v2begin, Iterator2 v2end);

/*! \relates Matrix */
template <class T> void swap(Matrix_t<T> &, Matrix_t<T> &);

/*! \relates Matrix */
template <class T>
std::ostream &operator<<(std::ostream &, const Matrix_t<T> &);

/*! \relates Matrix */
template <class T> Disposable<Matrix_t<T> > inverse(const Matrix_t<T> &m);

/*! \relates Matrix */
template <class T> T determinant(const Matrix_t<T> &m);

// inline definitions

template <class T>
inline Matrix_t<T>::Matrix_t()
    : data_((T *)(0)), rows_(0), columns_(0) {}

template <class T>
inline Matrix_t<T>::Matrix_t(Size rows, Size columns)
    : data_(rows * columns > 0 ? new T[rows * columns] : (T *)(0)), rows_(rows),
      columns_(columns) {}

template <class T>
inline Matrix_t<T>::Matrix_t(Size rows, Size columns, T value)
    : data_(rows * columns > 0 ? new T[rows * columns] : (T *)(0)), rows_(rows),
      columns_(columns) {
    std::fill(begin(), end(), value);
}

template<class T> template <class I>
Matrix_t<T>::Matrix_t(Size rows, Size columns, I start, I end)
    : data_(new T[rows * columns]), rows_(rows), columns_(columns) {
    std::copy(start, end, begin());
}


template <class T>
inline Matrix_t<T>::Matrix_t(const Matrix_t<T> &from)
    : data_(!from.empty() ? new T[from.rows_ * from.columns_] : (T *)(0)),
      rows_(from.rows_), columns_(from.columns_) {
#if defined(QL_PATCH_MSVC) && defined(QL_DEBUG)
    if (!from.empty())
#endif
        std::copy(from.begin(), from.end(), begin());
}

template <class T>
inline Matrix_t<T>::Matrix_t(const Disposable<Matrix_t<T> > &from)
    : data_((T *)(0)), rows_(0), columns_(0) {
    swap(const_cast<Disposable<Matrix_t<T> > &>(from));
}

template <class T>
inline Matrix_t<T> &Matrix_t<T>::operator=(const Matrix_t<T> &from) {
    // strong guarantee
    Matrix_t<T> temp(from);
    swap(temp);
    return *this;
}

template <class T>
inline Matrix_t<T> &Matrix_t<T>::
operator=(const Disposable<Matrix_t<T> > &from) {
    swap(const_cast<Disposable<Matrix_t<T> > &>(from));
    return *this;
}

template <class T> inline void Matrix_t<T>::swap(Matrix_t<T> &from) {
    using std::swap;
    data_.swap(from.data_);
    swap(rows_, from.rows_);
    swap(columns_, from.columns_);
}

template <class T>
inline const Matrix_t<T> &Matrix_t<T>::operator+=(const Matrix_t<T> &m) {
    QL_REQUIRE(rows_ == m.rows_ && columns_ == m.columns_,
               "matrices with different sizes (" << m.rows_ << "x" << m.columns_
                                                 << ", " << rows_ << "x"
                                                 << columns_ << ") cannot be "
                                                                "added");
    std::transform(begin(), end(), m.begin(), begin(), std::plus<T>());
    return *this;
}

template <class T>
inline const Matrix_t<T> &Matrix_t<T>::operator-=(const Matrix_t<T> &m) {
    QL_REQUIRE(rows_ == m.rows_ && columns_ == m.columns_,
               "matrices with different sizes (" << m.rows_ << "x" << m.columns_
                                                 << ", " << rows_ << "x"
                                                 << columns_ << ") cannot be "
                                                                "subtracted");
    std::transform(begin(), end(), m.begin(), begin(), std::minus<T>());
    return *this;
}

template <class T> inline const Matrix_t<T> &Matrix_t<T>::operator*=(T x) {
    std::transform(begin(), end(), begin(),
                   std::bind2nd(std::multiplies<T>(), x));
    return *this;
}

template <class T> inline const Matrix_t<T> &Matrix_t<T>::operator/=(T x) {
    std::transform(begin(), end(), begin(), std::bind2nd(std::divides<T>(), x));
    return *this;
}

template <class T>
inline typename Matrix_t<T>::const_iterator Matrix_t<T>::begin() const {
    return data_.get();
}

template <class T> inline typename Matrix_t<T>::iterator Matrix_t<T>::begin() {
    return data_.get();
}

template <class T>
inline typename Matrix_t<T>::const_iterator Matrix_t<T>::end() const {
    return data_.get() + rows_ * columns_;
}

template <class T> inline typename Matrix_t<T>::iterator Matrix_t<T>::end() {
    return data_.get() + rows_ * columns_;
}

template <class T>
inline typename Matrix_t<T>::const_reverse_iterator
Matrix_t<T>::rbegin() const {
    return const_reverse_iterator(end());
}

template <class T>
inline typename Matrix_t<T>::reverse_iterator Matrix_t<T>::rbegin() {
    return reverse_iterator(end());
}

template <class T>
inline typename Matrix_t<T>::const_reverse_iterator Matrix_t<T>::rend() const {
    return const_reverse_iterator(begin());
}

template <class T>
inline typename Matrix_t<T>::reverse_iterator Matrix_t<T>::rend() {
    return reverse_iterator(begin());
}

template <class T>
inline typename Matrix_t<T>::const_row_iterator
Matrix_t<T>::row_begin(Size i) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < rows_, "row index ("
                              << i << ") must be less than " << rows_
                              << ": matrix cannot be accessed out of range");
#endif
    return data_.get() + columns_ * i;
}

template <class T>
inline typename Matrix_t<T>::row_iterator Matrix_t<T>::row_begin(Size i) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < rows_, "row index ("
                              << i << ") must be less than " << rows_
                              << ": matrix cannot be accessed out of range");
#endif
    return data_.get() + columns_ * i;
}

template <class T>
inline typename Matrix_t<T>::const_row_iterator
Matrix_t<T>::row_end(Size i) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < rows_, "row index ("
                              << i << ") must be less than " << rows_
                              << ": matrix cannot be accessed out of range");
#endif
    return data_.get() + columns_ * (i + 1);
}

template <class T>
inline typename Matrix_t<T>::row_iterator Matrix_t<T>::row_end(Size i) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < rows_, "row index ("
                              << i << ") must be less than " << rows_
                              << ": matrix cannot be accessed out of range");
#endif
    return data_.get() + columns_ * (i + 1);
}

template <class T>
inline typename Matrix_t<T>::const_reverse_row_iterator
Matrix_t<T>::row_rbegin(Size i) const {
    return const_reverse_row_iterator(row_end(i));
}

template <class T>
inline typename Matrix_t<T>::reverse_row_iterator
Matrix_t<T>::row_rbegin(Size i) {
    return reverse_row_iterator(row_end(i));
}

template <class T>
inline typename Matrix_t<T>::const_reverse_row_iterator
Matrix_t<T>::row_rend(Size i) const {
    return const_reverse_row_iterator(row_begin(i));
}

template <class T>
inline typename Matrix_t<T>::reverse_row_iterator
Matrix_t<T>::row_rend(Size i) {
    return reverse_row_iterator(row_begin(i));
}

template <class T>
inline typename Matrix_t<T>::const_column_iterator
Matrix_t<T>::column_begin(Size i) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < columns_, "column index ("
                                 << i << ") must be less than " << columns_
                                 << ": matrix cannot be accessed out of range");
#endif
    return const_column_iterator(data_.get() + i, columns_);
}

template <class T>
inline typename Matrix_t<T>::column_iterator Matrix_t<T>::column_begin(Size i) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < columns_, "column index ("
                                 << i << ") must be less than " << columns_
                                 << ": matrix cannot be accessed out of range");
#endif
    return column_iterator(data_.get() + i, columns_);
}

template <class T>
inline typename Matrix_t<T>::const_column_iterator
Matrix_t<T>::column_end(Size i) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < columns_, "column index ("
                                 << i << ") must be less than " << columns_
                                 << ": matrix cannot be accessed out of range");
#endif
    return const_column_iterator(data_.get() + i + rows_ * columns_, columns_);
}

template <class T>
inline typename Matrix_t<T>::column_iterator Matrix_t<T>::column_end(Size i) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
    QL_REQUIRE(i < columns_, "column index ("
                                 << i << ") must be less than " << columns_
                                 << ": matrix cannot be accessed out of range");
#endif
    return column_iterator(data_.get() + i + rows_ * columns_, columns_);
}

template <class T>
inline typename Matrix_t<T>::const_reverse_column_iterator
Matrix_t<T>::column_rbegin(Size i) const {
    return const_reverse_column_iterator(column_end(i));
}

template <class T>
inline typename Matrix_t<T>::reverse_column_iterator
Matrix_t<T>::column_rbegin(Size i) {
    return reverse_column_iterator(column_end(i));
}

template <class T>
inline typename Matrix_t<T>::const_reverse_column_iterator
Matrix_t<T>::column_rend(Size i) const {
    return const_reverse_column_iterator(column_begin(i));
}

template <class T>
inline typename Matrix_t<T>::reverse_column_iterator
Matrix_t<T>::column_rend(Size i) {
    return reverse_column_iterator(column_begin(i));
}

template <class T>
inline typename Matrix_t<T>::const_row_iterator Matrix_t<T>::
operator[](Size i) const {
    return row_begin(i);
}

template <class T>
inline typename Matrix_t<T>::const_row_iterator Matrix_t<T>::at(Size i) const {
    QL_REQUIRE(i < rows_, "matrix access out of range");
    return row_begin(i);
}

template <class T>
inline typename Matrix_t<T>::row_iterator Matrix_t<T>::operator[](Size i) {
    return row_begin(i);
}

template <class T> inline typename Matrix_t<T>::row_iterator Matrix_t<T>::at(Size i) {
    QL_REQUIRE(i < rows_, "matrix access out of range");
    return row_begin(i);
}

template <class T>
inline Disposable<Array_t<T> > Matrix_t<T>::diagonal(void) const {
    Size arraySize = std::min<Size>(rows(), columns());
    Array_t<T> tmp(arraySize);
    for (Size i = 0; i < arraySize; i++)
        tmp[i] = (*this)[i][i];
    return tmp;
}

template <class T> inline Size Matrix_t<T>::rows() const { return rows_; }

template <class T> inline Size Matrix_t<T>::columns() const { return columns_; }

template <class T> inline bool Matrix_t<T>::empty() const {
    return rows_ == 0 || columns_ == 0;
}

template <class T>
inline const Disposable<Matrix_t<T> > operator+(const Matrix_t<T> &m1,
                                                const Matrix_t<T> &m2) {
    QL_REQUIRE(m1.rows() == m2.rows() && m1.columns() == m2.columns(),
               "matrices with different sizes ("
                   << m1.rows() << "x" << m1.columns() << ", " << m2.rows()
                   << "x" << m2.columns() << ") cannot be "
                                             "added");
    Matrix_t<T> temp(m1.rows(), m1.columns());
    std::transform(m1.begin(), m1.end(), m2.begin(), temp.begin(),
                   std::plus<T>());
    return temp;
}

template <class T>
inline const Disposable<Matrix_t<T> > operator-(const Matrix_t<T> &m1,
                                                const Matrix_t<T> &m2) {
    QL_REQUIRE(m1.rows() == m2.rows() && m1.columns() == m2.columns(),
               "matrices with different sizes ("
                   << m1.rows() << "x" << m1.columns() << ", " << m2.rows()
                   << "x" << m2.columns() << ") cannot be "
                                             "subtracted");
    Matrix_t<T> temp(m1.rows(), m1.columns());
    std::transform(m1.begin(), m1.end(), m2.begin(), temp.begin(),
                   std::minus<T>());
    return temp;
}

template <class T>
inline const Disposable<Matrix_t<T> > operator*(const Matrix_t<T> &m, T x) {
    Matrix_t<T> temp(m.rows(), m.columns());
    std::transform(m.begin(), m.end(), temp.begin(),
                   std::bind2nd(std::multiplies<T>(), x));
    return temp;
}

template <class T>
inline const Disposable<Matrix_t<T> > operator*(T x, const Matrix_t<T> &m) {
    Matrix_t<T> temp(m.rows(), m.columns());
    std::transform(m.begin(), m.end(), temp.begin(),
                   std::bind2nd(std::multiplies<T>(), x));
    return temp;
}

template <class T>
inline const Disposable<Matrix_t<T> > operator/(const Matrix_t<T> &m, T x) {
    Matrix_t<T> temp(m.rows(), m.columns());
    std::transform(m.begin(), m.end(), temp.begin(),
                   std::bind2nd(std::divides<T>(), x));
    return temp;
}

template <class T>
inline const Disposable<Array_t<T> > operator*(const Array_t<T> &v,
                                               const Matrix_t<T> &m) {
    QL_REQUIRE(v.size() == m.rows(),
               "vectors and matrices with different sizes ("
                   << v.size() << ", " << m.rows() << "x" << m.columns()
                   << ") cannot be multiplied");
    Array_t<T> result(m.columns());
    for (Size i = 0; i < result.size(); i++)
        result[i] =
            std::inner_product(v.begin(), v.end(), m.column_begin(i), T(0.0));
    return result;
}

template <class T>
inline const Disposable<Array_t<T> > operator*(const Matrix_t<T> &m,
                                               const Array_t<T> &v) {
    QL_REQUIRE(v.size() == m.columns(),
               "vectors and matrices with different sizes ("
                   << v.size() << ", " << m.rows() << "x" << m.columns()
                   << ") cannot be multiplied");
    Array_t<T> result(m.rows());
    for (Size i = 0; i < result.size(); i++)
        result[i] = std::inner_product(v.begin(), v.end(), m.row_begin(i), T(0.0));
    return result;
}

template <class T>
inline const Disposable<Matrix_t<T> > operator*(const Matrix_t<T> &m1,
                                                const Matrix_t<T> &m2) {
    QL_REQUIRE(m1.columns() == m2.rows(),
               "matrices with different sizes ("
                   << m1.rows() << "x" << m1.columns() << ", " << m2.rows()
                   << "x" << m2.columns() << ") cannot be "
                                             "multiplied");
    Matrix_t<T> result(m1.rows(), m2.columns(), 0.0);
    for (Size i = 0; i < result.rows(); ++i) {
        for (Size k = 0; k < m1.columns(); ++k) {
            for (Size j = 0; j < result.columns(); ++j) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return result;
}

template <class T>
inline const Disposable<Matrix_t<T> > transpose(const Matrix_t<T> &m) {
    Matrix_t<T> result(m.columns(), m.rows());
#if defined(QL_PATCH_MSVC) && defined(QL_DEBUG)
    if (!m.empty())
#endif
        for (Size i = 0; i < m.rows(); i++)
            std::copy(m.row_begin(i), m.row_end(i), result.column_begin(i));
    return result;
}

template <class T>
inline const Disposable<Matrix_t<T> > outerProduct(const Array_t<T> &v1,
                                                   const Array_t<T> &v2) {
    return outerProduct<T>(v1.begin(), v1.end(), v2.begin(), v2.end());
}

template <class T, class Iterator1, class Iterator2>
inline const Disposable<Matrix_t<T> >
outerProduct(Iterator1 v1begin, Iterator1 v1end, Iterator2 v2begin,
             Iterator2 v2end) {

    Size size1 = std::distance(v1begin, v1end);
    QL_REQUIRE(size1 > 0, "null first vector");

    Size size2 = std::distance(v2begin, v2end);
    QL_REQUIRE(size2 > 0, "null second vector");

    Matrix_t<T> result(size1, size2);

    for (Size i = 0; v1begin != v1end; i++, v1begin++)
        std::transform(v2begin, v2end, result.row_begin(i),
                       std::bind1st(std::multiplies<T>(), *v1begin));

    return result;
}

template <class T> inline void swap(Matrix_t<T> &m1, Matrix_t<T> &m2) {
    m1.swap(m2);
}

template <class T>
inline std::ostream &operator<<(std::ostream &out, const Matrix_t<T> &m) {
    std::streamsize width = out.width();
    for (Size i = 0; i < m.rows(); i++) {
        out << "| ";
        for (Size j = 0; j < m.columns(); j++)
            out << std::setw(int(width)) << m[i][j] << " ";
        out << "|\n";
    }
    return out;
}

typedef Matrix_t<Real> Matrix;

// implementation

template <class T> Disposable<Matrix_t<T> > inverse(const Matrix_t<T> &m) {
#if !defined(QL_NO_UBLAS_SUPPORT)

    QL_REQUIRE(m.rows() == m.columns(), "matrix is not square");

    boost::numeric::ublas::matrix<T> a(m.rows(), m.columns());

    std::copy(m.begin(), m.end(), a.data().begin());

    boost::numeric::ublas::permutation_matrix<Size> pert(m.rows());

    // lu decomposition
    const Size singular = lu_factorize(a, pert);
    QL_REQUIRE(singular == 0, "singular matrix given");

    boost::numeric::ublas::matrix<T> inverse =
        boost::numeric::ublas::identity_matrix<T>(m.rows());

    // backsubstitution
    boost::numeric::ublas::lu_substitute(a, pert, inverse);

    Matrix_t<T> retVal(m.rows(), m.columns());
    std::copy(inverse.data().begin(), inverse.data().end(), retVal.begin());

    return retVal;

#else
    QL_FAIL("this version of gcc does not support "
            "the Boost uBLAS library");
#endif
}

template <class T> T determinant(const Matrix_t<T> &m) {
#if !defined(QL_NO_UBLAS_SUPPORT)
    QL_REQUIRE(m.rows() == m.columns(), "matrix is not square");

    boost::numeric::ublas::matrix<T> a(m.rows(), m.columns());
    std::copy(m.begin(), m.end(), a.data().begin());

    // lu decomposition
    boost::numeric::ublas::permutation_matrix<Size> pert(m.rows());
    /* const Size singular = */ lu_factorize(a, pert);

    T retVal = 1.0;

    for (Size i = 0; i < m.rows(); ++i) {
        if (pert[i] != i)
            retVal *= -a(i, i);
        else
            retVal *= a(i, i);
    }
    return retVal;

#else
    QL_FAIL("this version of gcc does not support "
            "the Boost uBLAS library");
#endif
}

} // namespace QuantLib

#endif
