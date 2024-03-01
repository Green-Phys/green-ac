/*
 * Copyright (c) 2023 University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the “Software”), to deal in the Software
 * without restriction, including without limitation the rights to use, copy, modify,
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef GREEN_AC_NEVANLINNA_H
#define GREEN_AC_NEVANLINNA_H

#include <green/ndarray/ndarray.h>

#include <Eigen/Dense>
#include <complex>

#include "gmp_float.h"

namespace green::ac::nevanlinna {

  class nevanlinna_solver {
  public:
    using real_t    = gmp_float;
    using complex_t = std::complex<real_t>;
    using matrix_t  = Eigen::Matrix<complex_t, 2, 2>;
    /**
     * Build Nevanlinna factorization for data defined on positive Matsubara mesh
     *
     * @param mesh - Matsubara frequency mesh
     * @param data - Matsubara frequency data
     */
    void build(const std::vector<std::complex<double>>& mesh, const std::vector<std::complex<double>>& data);

    /**
     * Evaluate Nevanlinna continuation on a complex-valued grid
     * @param grid - grid for continuation evaluation
     * @return Nevanlinna continuation on the specified grid
     */
    [[nodiscard]] std::vector<std::complex<double>> evaluate(const std::vector<std::complex<double>>& grid);

  private:
    std::vector<complex_t>                          _phis;
    std::vector<matrix_t>                           _abcds;
    std::vector<complex_t>                          _mesh;
    std::vector<complex_t>                          _grid;
    std::vector<matrix_t>                           _coeffs;

    [[nodiscard]] std::vector<complex_t>            mobius_trasformation(const std::vector<std::complex<double>>& data) const;

    [[nodiscard]] std::vector<std::complex<double>> evaluate_internal(const std::vector<std::complex<double>>& grid) const;
  };

  class nevanlinna {
  public:
     nevanlinna() = default;
    ~nevanlinna() = default;

    // Copy/Move construction
    nevanlinna(nevanlinna const&) = default;
    nevanlinna(nevanlinna&&)      = default;

    /// Copy/Move assignment
    nevanlinna& operator=(nevanlinna const&) = default;
    nevanlinna& operator=(nevanlinna&&)      = default;

    /**
     * \brief Initialize Nevanlinna solvers for a given input data and solve Schur interpolation problem
     *
     * \param mesh - imaginary frequency mesh
     * \param data - data in imaginary frequency domain
     */
    void solve(const ndarray::ndarray<std::complex<double>, 1>& mesh, const ndarray::ndarray<std::complex<double>, 2>& data);

    /**
     * \brief Evaluate analytical continuation of imagiary frequency data using Nevanlinna method
     *
     * \param grid real frequency grid
     * \return analytically conrtinued data on a chosen real frequency grid
     */
    [[nodiscard]] ndarray::ndarray<std::complex<double>, 2> evaluate(std::vector<std::complex<double>>& grid);

  private:
    std::vector<nevanlinna_solver> _solvers{};
  };

}  // namespace green::ac::nevanlinna
#endif  // GREEN_AC_NEVANLINNA_H
