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

#include "green/ac/nevanlinna.h"

#include "green/ac/except.h"

namespace green::ac::nevanlinna {

  std::vector<nevanlinna_solver::complex_t> nevanlinna_solver::mobius_trasformation(
      const std::vector<std::complex<double>>& data) const {
    std::vector<complex_t> mdata(data.size());
    std::complex           I{0.0, 1.0};
    std::transform(data.begin(), data.end(), mdata.begin(),
                   [I](const std::complex<double>& d) { return complex_t(-d - I) / complex_t(-d + I); });
    return mdata;
  }

  void nevanlinna_solver::build(const std::vector<std::complex<double>>& mesh, const std::vector<std::complex<double>>& data) {
    assert(mesh.size() == data.size());
    if (std::any_of(mesh.begin(), mesh.end(), [](const std::complex<double>& v) { return v.real() != 0.0 or v.imag() < 0; })) {
      throw ac_nevanlinna_error("Data should be defined on the positive Matsubara frequencies.");
    }
    size_t M = mesh.size();
    _phis.resize(M);
    _abcds.resize(M);
    _mesh.resize(M);
    auto mdata = mobius_trasformation(data);
    _phis[0]   = mdata[0];
    for (int k = 0; k < M; k++) {
      _abcds[k] = matrix_t::Identity(2, 2);
      _mesh[k]  = mesh[k];
    }
    matrix_t prod(2, 2);
    for (int j = 0; j < M - 1; j++) {
      for (int k = j; k < M; k++) {
        prod << (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j])), _phis[j],
            std::conj(_phis[j]) * (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j])), complex_t{1., 0.};
        _abcds[k] *= prod;
      }
      _phis[j + 1] = (-_abcds[j + 1](1, 1) * mdata[j + 1] + _abcds[j + 1](0, 1)) /
                     (_abcds[j + 1](1, 0) * mdata[j + 1] - _abcds[j + 1](0, 0));
    }
    return;
  }

  std::vector<std::complex<double>> nevanlinna_solver::evaluate(const std::vector<std::complex<double>>& grid) {
    size_t M = _phis.size();
    if (M == 0) {
      throw ac_nevanlinna_error("Empty continuation data. Please run solve(...) first.");
    }
    if (grid.size() == _grid.size() &&
        std::equal(grid.begin(), grid.end(), _grid.begin(), [](const std::complex<double>& w1, const complex_t& w2) {
          return std::abs(w1.real() - w2.real().convert_to<double>()) < 1e-9 &&
                 std::abs(w1.imag() - w2.imag().convert_to<double>()) < 1e-9;
        })) {
      return evaluate_internal(grid);
    }
    std::cout << "Reevaluation" << std::endl;
    _coeffs = std::vector<matrix_t>(grid.size());
    _grid.resize(grid.size());
    std::transform(grid.begin(), grid.end(), _grid.begin(), [](const std::complex<double>& w) {
      return complex_t{w.real(), w.imag()};
    });
    complex_t                         I{0., 1.};
    complex_t                         One{1., 0.};
    std::vector<std::complex<double>> results(grid.size());
    matrix_t                          prod(2, 2);
    for (int i = 0; i < grid.size(); ++i) {
      matrix_t result = matrix_t::Identity(2, 2);
      auto     z      = complex_t(grid[i]);
      for (int j = 0; j < M; j++) {
        prod << (z - _mesh[j]) / (z - std::conj(_mesh[j])), _phis[j],
            std::conj(_phis[j]) * ((z - _mesh[j]) / (z - std::conj(_mesh[j]))), complex_t{1., 0.};
        result *= prod;
      }
      _coeffs[i] = result;
      complex_t param{0., 0.};
      complex_t theta = (result(0, 0) * param + result(0, 1)) / (result(1, 0) * param + result(1, 1));
      complex_t value = I * (One + theta) / (One - theta);
      results[i]      = -std::complex<double>(value.real().convert_to<double>(),
                                         value.imag().convert_to<double>());  // inverse Mobius transform from theta to NG
    }
    return results;
  }

  std::vector<std::complex<double>> nevanlinna_solver::evaluate_internal(const std::vector<std::complex<double>>& grid) const {
    size_t M = _phis.size();
    if (M == 0) {
      throw ac_nevanlinna_error("Empty continuation data. Please run solve(...) first.");
    }
    complex_t                         I{0., 1.};
    complex_t                         One{1., 0.};
    std::vector<std::complex<double>> results(grid.size());
    for (int i = 0; i < _grid.size(); ++i) {
      matrix_t  result         = _coeffs[i];
      auto      z              = _grid[i];
      complex_t theta_M_plus_1 = complex_t{0., 0.};

      complex_t theta          = (result(0, 0) * theta_M_plus_1 + result(0, 1)) / (result(1, 0) * theta_M_plus_1 + result(1, 1));
      complex_t value          = I * (One + theta) / (One - theta);
      // inverse Mobius transform from theta to NG
      results[i] = -std::complex<double>(value.real().convert_to<double>(), value.imag().convert_to<double>());
    }
    return results;
  }

  void nevanlinna::solve(const ndarray::ndarray<std::complex<double>, 1>& mesh,
                         const ndarray::ndarray<std::complex<double>, 2>& data) {
    _solvers.clear();
    size_t N  = data.shape()[1];
    size_t nw = std::count_if(mesh.begin(), mesh.end(), [](const std::complex<double>& w) { return w.imag() > 0; });
    std::vector<std::complex<double>> data_in(nw);
    std::vector<std::complex<double>> mesh_in(nw);
    for (size_t n = 0; n < N; ++n) {
      for (size_t iw = 0, iww = 0; iw < mesh.shape()[0]; ++iw) {
        if (mesh(iw).imag() < 0) continue;
        data_in[iww] = data(iw, n);
        mesh_in[iww] = mesh(iw);
        ++iww;
      }
      nevanlinna_solver f;
      f.build(mesh_in, data_in);
      _solvers.push_back(f);
    }
  }

  ndarray::ndarray<std::complex<double>, 2> nevanlinna::evaluate(std::vector<std::complex<double>>& grid) {
    size_t                                    size = _solvers.size();
    ndarray::ndarray<std::complex<double>, 2> results(grid.size(), size);
    for (size_t n = 0; n < size; ++n) {
      std::vector<std::complex<double>> data = _solvers[n].evaluate(grid);
      for (size_t iw = 0; iw < grid.size(); ++iw) {
        results(iw, n) = data[iw];
      }
    }
    return results;
  }
}  // namespace green::ac::nevanlinna
