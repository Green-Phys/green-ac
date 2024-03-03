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

#include <green/ac/except.h>
#include <green/ac/nevanlinna.h>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

using namespace green::ac;
using namespace green::ndarray;

TEST_CASE("Nevanlinna") {
  const double                     eta       = 0.1;
  const double                     mu1       = 1.0;
  const double                     beta      = 10.0;
  const double                     omega_min = -5.0;
  const double                     omega_max = 5.0;
  const int                        n_iw      = 100;
  const int                        n_omega   = 100;
  nevanlinna::nevanlinna           a;
  ndarray<std::complex<double>, 1> im_data(n_iw);
  ndarray<std::complex<double>, 1> im_grid(n_iw);
  ndarray<std::complex<double>, 1> re_data(n_omega);
  ndarray<std::complex<double>, 1> re_grid(n_omega);

  for (int iw = 0, w = -n_iw / 2; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2 * w + 1) * M_PI * std::complex<double>(0, 1) / beta;
    im_data(iw) = 1. / (im_grid(iw) - mu1);
  }

  for (size_t iw = 0; iw < n_omega; ++iw) {
    re_grid(iw) = std::complex<double>(iw * (omega_max - omega_min) / n_omega, eta);
    re_data(iw) = 1. / (re_grid(iw) - mu1);
  }

  a.solve(im_grid, im_data);
  ndarray<std::complex<double>, 1> result = a.evaluate(re_grid);

  REQUIRE(std::equal(result.begin(), result.end(), re_data.begin(),
                     [](const std::complex<double>& x, const std::complex<double>& y) { return std::abs(x - y) < 1e-12; }));
}

TEST_CASE("Nevanlinna Solver") {
  const double                     eta       = 0.1;
  const double                     mu        = 1.0;
  const double                     beta      = 10.0;
  const double                     omega_min = -5.0;
  const double                     omega_max = 5.0;
  const int                        n_iw      = 100;
  const int                        n_omega   = 100;
  nevanlinna::nevanlinna           a;
  ndarray<std::complex<double>, 1> im_data(n_iw);
  ndarray<std::complex<double>, 1> im_grid(n_iw);
  ndarray<std::complex<double>, 1> re_data(n_omega);
  ndarray<std::complex<double>, 1> re_grid(n_omega);

  for (int iw = 0, w = -n_iw / 2; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2 * w + 1) * M_PI * std::complex<double>(0, 1) / beta;
    im_data(iw) = 1. / (im_grid(iw) - mu);
  }

  REQUIRE_THROWS_AS(a.build(im_grid, im_data), ac_nevanlinna_error);
  REQUIRE_THROWS_AS(a.evaluate(re_grid), ac_nevanlinna_error);

  for (int iw = 0, w = -n_iw / 2; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2 * iw + 1) * M_PI * std::complex<double>(0, 1) / beta;
    im_data(iw) = 1. / (im_grid(iw) - mu);
  }

  for (size_t iw = 0; iw < n_omega; ++iw) {
    re_grid(iw) = std::complex<double>(iw * (omega_max - omega_min) / n_omega, eta);
    re_data(iw) = 1. / (re_grid(iw) - mu);
  }

  a.build(im_grid, im_data);
  ndarray<std::complex<double>, 1> result = a.evaluate(re_grid);

  REQUIRE(std::equal(result.begin(), result.end(), re_data.begin(),
                     [](const std::complex<double>& x, const std::complex<double>& y) { return std::abs(x - y) < 1e-12; }));
  ndarray<std::complex<double>, 1> result2 = a.evaluate(re_grid);
  REQUIRE(std::equal(result2.begin(), result2.end(), re_data.begin(),
                     [](const std::complex<double>& x, const std::complex<double>& y) { return std::abs(x - y) < 1e-12; }));
}