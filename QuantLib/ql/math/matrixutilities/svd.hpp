/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Neil Firth
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

    Adapted from the TNT project
    http://math.nist.gov/tnt/download.html

    This software was developed at the National Institute of Standards
    and Technology (NIST) by employees of the Federal Government in the
    course of their official duties. Pursuant to title 17 Section 105
    of the United States Code this software is not subject to copyright
    protection and is in the public domain. NIST assumes no responsibility
    whatsoever for its use by other parties, and makes no guarantees,
    expressed or implied, about its quality, reliability, or any other
    characteristic.

    We would appreciate acknowledgement if the software is incorporated in
    redistributable libraries or applications.

*/

/*! \file svd.hpp
    \brief singular value decomposition
*/

#ifndef quantlib_math_svd_h
#define quantlib_math_svd_h

#include <ql/math/matrix.hpp>

namespace QuantLib {

//! Singular value decomposition
/*! Refer to Golub and Van Loan: Matrix computation,
    The Johns Hopkins University Press

    \test the correctness of the returned values is tested by
          checking their properties.
*/
template <class T> class SVD_t {
  public:
    // constructor
    SVD_t(const Matrix_t<T> &);
    // results
    const Matrix_t<T> &U() const;
    const Matrix_t<T> &V() const;
    const Array_t<T> &singularValues() const;
    Disposable<Matrix_t<T> > S() const;
    T norm2() const;
    T cond() const;
    Size rank() const;
    // utilities
    Disposable<Array_t<T> > solveFor(const Array_t<T> &) const;

  private:
    Matrix_t<T> U_, V_;
    Array_t<T> s_;
    Integer m_, n_;
    bool transpose_;
};

typedef SVD_t<Real> SVD;

// implementation

namespace {

/*  returns hypotenuse of real (non-complex) scalars a and b by
    avoiding underflow/overflow
    using (a * sqrt( 1 + (b/a) * (b/a))), rather than
    sqrt(a*a + b*b).
*/
template <class T> inline T hypot(const T &a, const T &b) {
    if (a == 0) {
        return QLFCT::abs(b);
    } else {
        T c = b / a;
        return QLFCT::abs(a) * QLFCT::sqrt(1 + c * c);
    }
}
}

template <class T> SVD_t<T>::SVD_t(const Matrix_t<T> &M) {

    using std::swap;

    Matrix_t<T> A;

    /* The implementation requires that rows > columns.
       If this is not the case, we decompose M^T instead.
       Swapping the resulting U and V gives the desired
       result for M as

       M^T = U S V^T           (decomposition of M^T)

       M = (U S V^T)^T         (transpose)

       M = (V^T^T S^T U^T)     ((AB)^T = B^T A^T)

       M = V S U^T             (idempotence of transposition,
                                symmetry of diagonal matrix S)

    */

    if (M.rows() >= M.columns()) {
        A = M;
        transpose_ = false;
    } else {
        A = transpose(M);
        transpose_ = true;
    }

    m_ = A.rows();
    n_ = A.columns();

    // we're sure that m_ >= n_

    s_ = Array_t<T>(n_);
    U_ = Matrix_t<T>(m_, n_, 0.0);
    V_ = Matrix_t<T>(n_, n_);
    Array_t<T> e(n_);
    Array_t<T> work(m_);
    Integer i, j, k;

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.

    Integer nct = static_cast<Integer>(QLFCT::min<double>(m_ - 1, n_));
    Integer nrt = static_cast<Integer>(QLFCT::max<double>(0, n_ - 2));
    for (k = 0; k < static_cast<Integer>(QLFCT::max<double>(nct, nrt)); k++) {
        if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            s_[k] = 0;
            for (i = k; i < m_; i++) {
                s_[k] = hypot(s_[k], A[i][k]);
            }
            if (s_[k] != 0.0) {
                if (A[k][k] < 0.0) {
                    s_[k] = -s_[k];
                }
                for (i = k; i < m_; i++) {
                    A[i][k] /= s_[k];
                }
                A[k][k] += 1.0;
            }
            s_[k] = -s_[k];
        }
        for (j = k + 1; j < n_; j++) {
            if ((k < nct) && (s_[k] != 0.0)) {

                // Apply the transformation.

                T t = 0;
                for (i = k; i < m_; i++) {
                    t += A[i][k] * A[i][j];
                }
                t = -t / A[k][k];
                for (i = k; i < m_; i++) {
                    A[i][j] += t * A[i][k];
                }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e[j] = A[k][j];
        }
        if (k < nct) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (i = k; i < m_; i++) {
                U_[i][k] = A[i][k];
            }
        }
        if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (i = k + 1; i < n_; i++) {
                e[k] = hypot(e[k], e[i]);
            }
            if (e[k] != 0.0) {
                if (e[k + 1] < 0.0) {
                    e[k] = -e[k];
                }
                for (i = k + 1; i < n_; i++) {
                    e[i] /= e[k];
                }
                e[k + 1] += 1.0;
            }
            e[k] = -e[k];
            if ((k + 1 < m_) & (e[k] != 0.0)) {

                // Apply the transformation.

                for (i = k + 1; i < m_; i++) {
                    work[i] = 0.0;
                }
                for (j = k + 1; j < n_; j++) {
                    for (i = k + 1; i < m_; i++) {
                        work[i] += e[j] * A[i][j];
                    }
                }
                for (j = k + 1; j < n_; j++) {
                    T t = -e[j] / e[k + 1];
                    for (i = k + 1; i < m_; i++) {
                        A[i][j] += t * work[i];
                    }
                }
            }

            // Place the transformation in V for subsequent
            // back multiplication.

            for (i = k + 1; i < n_; i++) {
                V_[i][k] = e[i];
            }
        }
    }

    // Set up the final bidiagonal matrix or order n.

    if (nct < n_) {
        s_[nct] = A[nct][nct];
    }
    if (nrt + 1 < n_) {
        e[nrt] = A[nrt][n_ - 1];
    }
    e[n_ - 1] = 0.0;

    // generate U

    for (j = nct; j < n_; j++) {
        for (i = 0; i < m_; i++) {
            U_[i][j] = 0.0;
        }
        U_[j][j] = 1.0;
    }
    for (k = nct - 1; k >= 0; --k) {
        if (s_[k] != 0.0) {
            for (j = k + 1; j < n_; ++j) {
                T t = 0;
                for (i = k; i < m_; i++) {
                    t += U_[i][k] * U_[i][j];
                }
                t = -t / U_[k][k];
                for (i = k; i < m_; i++) {
                    U_[i][j] += t * U_[i][k];
                }
            }
            for (i = k; i < m_; i++) {
                U_[i][k] = -U_[i][k];
            }
            U_[k][k] = 1.0 + U_[k][k];
            for (i = 0; i < k - 1; i++) {
                U_[i][k] = 0.0;
            }
        } else {
            for (i = 0; i < m_; i++) {
                U_[i][k] = 0.0;
            }
            U_[k][k] = 1.0;
        }
    }

    // generate V

    for (k = n_ - 1; k >= 0; --k) {
        if ((k < nrt) & (e[k] != 0.0)) {
            for (j = k + 1; j < n_; ++j) {
                T t = 0;
                for (i = k + 1; i < n_; i++) {
                    t += V_[i][k] * V_[i][j];
                }
                t = -t / V_[k + 1][k];
                for (i = k + 1; i < n_; i++) {
                    V_[i][j] += t * V_[i][k];
                }
            }
        }
        for (i = 0; i < n_; i++) {
            V_[i][k] = 0.0;
        }
        V_[k][k] = 1.0;
    }

    // Main iteration loop for the singular values.

    Integer p = n_, pp = p - 1;
    Integer iter = 0;
    T eps = QLFCT::pow(2.0, -52.0);
    while (p > 0) {
        Integer k;
        Integer kase;

        // Here is where a test for too many iterations would go.

        // This section of the program inspects for
        // negligible elements in the s and e arrays.  On
        // completion the variables kase and k are set as follows.

        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e[k-1] is negligible, k<p, and
        //              s(k), ..., s(p) are not negligible (qr step).
        // kase = 4     if e(p-1) is negligible (convergence).

        for (k = p - 2; k >= -1; --k) {
            if (k == -1) {
                break;
            }
            if (QLFCT::abs(e[k]) <=
                eps * (QLFCT::abs(s_[k]) + QLFCT::abs(s_[k + 1]))) {
                e[k] = 0.0;
                break;
            }
        }
        if (k == p - 2) {
            kase = 4;
        } else {
            Integer ks;
            for (ks = p - 1; ks >= k; --ks) {
                if (ks == k) {
                    break;
                }
                T t = (ks != p ? QLFCT::abs(e[ks]) : 0.) +
                      (ks != k + 1 ? QLFCT::abs(e[ks - 1]) : 0.);
                if (QLFCT::abs(s_[ks]) <= eps * t) {
                    s_[ks] = 0.0;
                    break;
                }
            }
            if (ks == k) {
                kase = 3;
            } else if (ks == p - 1) {
                kase = 1;
            } else {
                kase = 2;
                k = ks;
            }
        }
        k++;

        // Perform the task indicated by kase.

        switch (kase) {

        // Deflate negligible s(p).

        case 1: {
            T f = e[p - 2];
            e[p - 2] = 0.0;
            for (j = p - 2; j >= k; --j) {
                T t = hypot(s_[j], f);
                T cs = s_[j] / t;
                T sn = f / t;
                s_[j] = t;
                if (j != k) {
                    f = -sn * e[j - 1];
                    e[j - 1] = cs * e[j - 1];
                }
                for (i = 0; i < n_; i++) {
                    t = cs * V_[i][j] + sn * V_[i][p - 1];
                    V_[i][p - 1] = -sn * V_[i][j] + cs * V_[i][p - 1];
                    V_[i][j] = t;
                }
            }
        } break;

        // Split at negligible s(k).

        case 2: {
            T f = e[k - 1];
            e[k - 1] = 0.0;
            for (j = k; j < p; j++) {
                T t = hypot(s_[j], f);
                T cs = s_[j] / t;
                T sn = f / t;
                s_[j] = t;
                f = -sn * e[j];
                e[j] = cs * e[j];
                for (i = 0; i < m_; i++) {
                    t = cs * U_[i][j] + sn * U_[i][k - 1];
                    U_[i][k - 1] = -sn * U_[i][j] + cs * U_[i][k - 1];
                    U_[i][j] = t;
                }
            }
        } break;

        // Perform one qr step.

        case 3: {

            // Calculate the shift.
            T scale = QLFCT::max(
                QLFCT::max(QLFCT::max(QLFCT::max(QLFCT::abs(s_[p - 1]),
                                                 QLFCT::abs(s_[p - 2])),
                                      QLFCT::abs(e[p - 2])),
                           QLFCT::abs(s_[k])),
                QLFCT::abs(e[k]));
            T sp = s_[p - 1] / scale;
            T spm1 = s_[p - 2] / scale;
            T epm1 = e[p - 2] / scale;
            T sk = s_[k] / scale;
            T ek = e[k] / scale;
            T b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
            T c = (sp * epm1) * (sp * epm1);
            T shift = 0.0;
            if ((b != 0.0) | (c != 0.0)) {
                shift = QLFCT::sqrt(b * b + c);
                if (b < 0.0) {
                    shift = -shift;
                }
                shift = c / (b + shift);
            }
            T f = (sk + sp) * (sk - sp) + shift;
            T g = sk * ek;

            // Chase zeros.

            for (j = k; j < p - 1; j++) {
                T t = hypot(f, g);
                T cs = f / t;
                T sn = g / t;
                if (j != k) {
                    e[j - 1] = t;
                }
                f = cs * s_[j] + sn * e[j];
                e[j] = cs * e[j] - sn * s_[j];
                g = sn * s_[j + 1];
                s_[j + 1] = cs * s_[j + 1];
                for (i = 0; i < n_; i++) {
                    t = cs * V_[i][j] + sn * V_[i][j + 1];
                    V_[i][j + 1] = -sn * V_[i][j] + cs * V_[i][j + 1];
                    V_[i][j] = t;
                }
                t = hypot(f, g);
                cs = f / t;
                sn = g / t;
                s_[j] = t;
                f = cs * e[j] + sn * s_[j + 1];
                s_[j + 1] = -sn * e[j] + cs * s_[j + 1];
                g = sn * e[j + 1];
                e[j + 1] = cs * e[j + 1];
                if (j < m_ - 1) {
                    for (i = 0; i < m_; i++) {
                        t = cs * U_[i][j] + sn * U_[i][j + 1];
                        U_[i][j + 1] = -sn * U_[i][j] + cs * U_[i][j + 1];
                        U_[i][j] = t;
                    }
                }
            }
            e[p - 2] = f;
            iter = iter + 1;
        } break;

        // Convergence.

        case 4: {

            // Make the singular values positive.

            if (s_[k] <= 0.0) {
                s_[k] = (s_[k] < 0.0 ? -s_[k] : 0.0);
                for (i = 0; i <= pp; i++) {
                    V_[i][k] = -V_[i][k];
                }
            }

            // Order the singular values.

            while (k < pp) {
                if (s_[k] >= s_[k + 1]) {
                    break;
                }
                swap(s_[k], s_[k + 1]);
                if (k < n_ - 1) {
                    for (i = 0; i < n_; i++) {
                        swap(V_[i][k], V_[i][k + 1]);
                    }
                }
                if (k < m_ - 1) {
                    for (i = 0; i < m_; i++) {
                        swap(U_[i][k], U_[i][k + 1]);
                    }
                }
                k++;
            }
            iter = 0;
            --p;
        } break;
        }
    }
}

template <class T> const Matrix_t<T> &SVD_t<T>::U() const {
    return (transpose_ ? V_ : U_);
}

template <class T> const Matrix_t<T> &SVD_t<T>::V() const {
    return (transpose_ ? U_ : V_);
}

template <class T> const Array_t<T> &SVD_t<T>::singularValues() const {
    return s_;
}

template <class T> Disposable<Matrix_t<T> > SVD_t<T>::S() const {
    Matrix_t<T> S(n_, n_);
    for (Size i = 0; i < Size(n_); i++) {
        for (Size j = 0; j < Size(n_); j++) {
            S[i][j] = 0.0;
        }
        S[i][i] = s_[i];
    }
    return S;
}

template <class T> T SVD_t<T>::norm2() const { return s_[0]; }

template <class T> T SVD_t<T>::cond() const { return s_[0] / s_[n_ - 1]; }

template <class T> Size SVD_t<T>::rank() const {
    T eps = QL_EPSILON;
    T tol = m_ * s_[0] * eps;
    Size r = 0;
    for (Size i = 0; i < s_.size(); i++) {
        if (s_[i] > tol) {
            r++;
        }
    }
    return r;
}

template <class T>
Disposable<Array_t<T> > SVD_t<T>::solveFor(const Array_t<T> &b) const {
    Matrix_t<T> W(n_, n_, 0.0);
    const Size numericalRank = this->rank();
    for (Size i = 0; i < numericalRank; i++)
        W[i][i] = 1. / s_[i];

    Matrix_t<T> inverse = V() * W * transpose(U());
    Array_t<T> result = inverse * b;
    return result;
}
}

#endif
