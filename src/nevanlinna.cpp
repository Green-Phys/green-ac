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

#include <fstream>

#include "green/ac/except.h"
namespace green::ac::nevanlinna {

  std::vector<nevanlinna::complex_t> nevanlinna::mobius_trasformation(const array_t& data) const {
    std::vector<complex_t> mdata(data.size());
    std::complex           I{0.0, 1.0};
    std::transform(data.begin(), data.end(), mdata.begin(),
                   [I](const std::complex<double>& d) { return complex_t(-d - I) / complex_t(-d + I); });
    return mdata;
  }  // LCOV_EXCL_LINE

  void nevanlinna::build(const array_t& mesh, const array_t& data) {
    assert(mesh.size() == data.size());
    mpf_set_default_prec(_precision);
    if (std::any_of(mesh.begin(), mesh.end(), [](const std::complex<double>& v) { return v.real() != 0.0 or v.imag() < 0; })) {
      throw ac_nevanlinna_error("Data should be defined on the positive Matsubara frequencies.");
    }
    size_t M = mesh.size();
    _phis.resize(M);
    _abcds.resize(M);
    _mesh.resize(M);
    auto mdata = mobius_trasformation(data);
    for (int k = 0; k < M; k++) {
      _abcds[k] = matrix_t::Identity(2, 2);
      _mesh[k]  = mesh(k);
    }
    std::reverse(mdata.begin(), mdata.end());
    std::reverse(_mesh.begin(), _mesh.end());
    _phis[0]   = mdata[0];
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
  }

  nevanlinna::array_t nevanlinna::evaluate(const array_t& grid) {
    mpf_set_default_prec(_precision);
    size_t M = _phis.size();
    if (M == 0) {
      throw ac_nevanlinna_error("Empty continuation data. Please run solve(...) first.");
    }
    if (grid.size() == _grid.size() &&
        std::equal(grid.begin(), grid.end(), _grid.begin(), [](const std::complex<double>& w1, const complex_t& w2) {
          return std::abs(w1.real() - to_double(w2.real())) < 1e-9 && std::abs(w1.imag() - to_double(w2.imag())) < 1e-9;
        })) {
      return evaluate_internal(grid);
    }
    _coeffs = std::vector<matrix_t>(grid.size());
    _grid.resize(grid.size());
    std::transform(grid.begin(), grid.end(), _grid.begin(), [](const std::complex<double>& w) {
      return complex_t{w.real(), w.imag()};
    });
    complex_t     I{0., 1.};
    complex_t     One{1., 0.};
    array_t       results(grid.size());
    matrix_t      prod(2, 2);
    for (int i = 0; i < grid.size(); ++i) {
      matrix_t result = matrix_t::Identity(2, 2);
      auto     z      = _grid[i];
      for (int j = 0; j < M; j++) {
        prod << (z - _mesh[j]) / (z - std::conj(_mesh[j])), _phis[j],
            std::conj(_phis[j]) * ((z - _mesh[j]) / (z - std::conj(_mesh[j]))), One;
        result *= prod;
      }
      _coeffs[i] = result;
      complex_t param{0., 0.};
      complex_t theta = (result(0, 0) * param + result(0, 1)) / (result(1, 0) * param + result(1, 1));
      // inverse Mobius transform from theta to NG
      complex_t value = I * (One + theta) / (One - theta);
      results(i)      = -std::complex<double>(to_double(value.real()), to_double(value.imag()));
    }
    return results;
  }

  nevanlinna::array_t nevanlinna::evaluate_internal(const array_t& grid) const {
    complex_t I{0., 1.};
    complex_t One{1., 0.};
    array_t   results(grid.size());
    for (int i = 0; i < _grid.size(); ++i) {
      matrix_t  result         = _coeffs[i];
      auto      z              = _grid[i];
      complex_t theta_M_plus_1 = complex_t{0., 0.};

      complex_t theta          = (result(0, 0) * theta_M_plus_1 + result(0, 1)) / (result(1, 0) * theta_M_plus_1 + result(1, 1));
      complex_t value          = I * (One + theta) / (One - theta);
      // inverse Mobius transform from theta to NG
      results(i) = -std::complex<double>(to_double(value.real()), to_double(value.imag()));
    }
    return results;
  }

  void nevanlinna::solve(const array_t& mesh, const array_t& data) {
    size_t  nw = std::count_if(mesh.begin(), mesh.end(), [](const std::complex<double>& w) { return w.imag() > 0; });
    // invalidate solved data
    _grid.resize(0);
    array_t data_in(nw);
    array_t mesh_in(nw);
    for (int64_t iw = 0, iww = 0; iw < mesh.shape()[0]; ++iw) {
      if (mesh(iw).imag() < 0) continue;
      data_in(iww) = data(iw);
      mesh_in(iww) = mesh(iw);
      ++iww;
    }
    build(mesh_in, data_in);
  }
}  // namespace green::ac::nevanlinna
