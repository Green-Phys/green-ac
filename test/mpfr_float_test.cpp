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

#include <green/ac/mpfr_float.h>

#include <Eigen/Dense>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

using namespace green::ac;
TEST_CASE("MPFR Math") {
  SECTION("Inplace") {
    mpfr_float a(10);
    mpfr_float b(20);
    a += b;
    REQUIRE(std::abs(double(a) - 30.0) < 1e-12);
    a -= mpfr_float(30);
    REQUIRE(std::abs(double(a)) < 1e-12);
    b /= 4.0;
    REQUIRE(std::abs(double(b) - 5.0) < 1e-12);
    b *= mpfr_float(2);
    REQUIRE(std::abs(double(b) - 10.0) < 1e-12);
  }

  SECTION("Normal") {
    mpfr_float a(10);
    mpfr_float b(20);
    REQUIRE(std::abs(double(a + b) - 30.0) < 1e-12);
    REQUIRE(std::abs(double(a - b) - -10.0) < 1e-12);
    REQUIRE(std::abs(double(a * b) - 200.0) < 1e-12);
    REQUIRE(std::abs(double(b / a) - 2.0) < 1e-12);
    REQUIRE(std::abs(double(cos(a)) - std::cos(double(a))) < 1e-12);
    REQUIRE(std::abs(double(sin(a)) - std::sin(double(a))) < 1e-12);
    REQUIRE(std::abs(double(sqrt(a)) - std::sqrt(double(a))) < 1e-12);
    REQUIRE(std::abs(double(atan2(a, b)) - std::atan2(double(a), double(b))) < 1e-12);
  }

  SECTION("Mixed types") {
    mpfr_float a(10);
    REQUIRE(std::abs(double(a + 10) - 20.0) < 1e-12);
    REQUIRE(std::abs(double(a - 20) - -10.0) < 1e-12);
    REQUIRE(std::abs(double(a * 5) - 50.0) < 1e-12);
    REQUIRE(std::abs(double(a / 4) - 2.5) < 1e-12);
    REQUIRE(std::abs(double(10 + a) - 20.0) < 1e-12);
    REQUIRE(std::abs(double(20 - a) - 10.0) < 1e-12);
    REQUIRE(std::abs(double(5 * a) - 50.0) < 1e-12);
    REQUIRE(std::abs(double(4 / a) - 0.4) < 1e-12);
  }
}

TEST_CASE("MPFR Complex") {
  SECTION("Init") {
    std::complex             a(mpfr_float(1), mpfr_float(1));
    std::complex<mpfr_float> b(1.0);
    std::complex<mpfr_float> c(std::complex(1.0, 1.0));
    REQUIRE(std::abs(double(a.imag()) - 1.0) < 1e-12);
    REQUIRE(std::abs(double(b.real()) - 1.0) < 1e-12);
    REQUIRE(std::abs(double(b.imag())) < 1e-12);
    REQUIRE(std::abs(double(c.real()) - 1.0) < 1e-12);
    REQUIRE(std::abs(double(c.imag()) - 1.0) < 1e-12);
    a = b;
    REQUIRE(std::abs(double(a.imag())) < 1e-12);
    c = 5;
    REQUIRE(std::abs(double(c.real() - 5)) < 1e-12);
    REQUIRE(std::abs(double(c.imag())) < 1e-12);
    a += std::complex{2., 3.};
    REQUIRE(std::abs(double(a.real() - 3.0)) < 1e-12);
    REQUIRE(std::abs(double(a.imag() - 3.0)) < 1e-12);
    a -= std::complex{2., 3.};
    REQUIRE(std::abs(double(a.real() - 1.0)) < 1e-12);
    REQUIRE(std::abs(double(a.imag())) < 1e-12);
  }
  SECTION("Math") {
    std::complex a(mpfr_float(1), mpfr_float(1));
    std::complex b(mpfr_float(6), mpfr_float(3));
    std::complex ad(1, 1);
    std::complex bd(6, 3);

    auto         c = a * b;
    REQUIRE(std::abs(double(c.real()) - 3.0) < 1e-12);
    REQUIRE(std::abs(double(c.imag()) - 9.0) < 1e-12);

    REQUIRE(std::abs(double(std::abs(std::complex{mpfr_float(2.), mpfr_float(2.)}) - std::abs(std::complex{2., 2.}))) < 1e-12);
    REQUIRE(double(std::abs(std::polar(mpfr_float(2.), mpfr_float(2.)) - std::polar(2., 2.))) < 1e-12);
    REQUIRE(double(std::abs(std::conj(std::complex{mpfr_float(2.), mpfr_float(2.)}) - std::conj(std::complex{2., 2.}))) < 1e-12);
    REQUIRE(std::abs(double(
                std::abs(std::sqrt(std::complex{mpfr_float(2.), mpfr_float(2.)}) - std::sqrt(std::complex{2., 2.})))) < 1e-12);
  }
  SECTION("Mixed Type Math") {
    std::complex a(mpfr_float(1), mpfr_float(1));
    std::complex b(mpfr_float(6), mpfr_float(3));
    std::complex ad(1., 1.);
    std::complex bd(6., 3.);
    a += 1.0;
    REQUIRE(std::abs(double(a.real()) - 2.0) < 1e-12);
    a -= 1;
    REQUIRE(std::abs(double(a.real()) - 1.0) < 1e-12);
    auto c = a + std::complex(0.0, 1.0);
    REQUIRE(std::abs(double(c.imag()) - 2.0) < 1e-12);
    auto a_ad = a * ad;
    REQUIRE(std::abs(double(a_ad.real())) < 1e-12);
    REQUIRE(std::abs(double(a_ad.imag()) - 2.0) < 1e-12);
    REQUIRE(std::abs(double(std::abs(a / b)) - std::abs(ad / bd)) < 1e-12);
  }

  SECTION("Matrix") {
    using MatrixTf = Eigen::Matrix<mpfr_float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MatrixTz = Eigen::Matrix<std::complex<mpfr_float>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    MatrixTf a(2, 2);
    MatrixTf b(2, 2);
    a << mpfr_float(2), mpfr_float(0), mpfr_float(0), mpfr_float(2);
    b << mpfr_float(2), mpfr_float(0), mpfr_float(0), mpfr_float(2);
    auto c = a * b;
    REQUIRE(double(c(0, 0)) == 4.0);
    a *= b;
    REQUIRE(double(a(0, 0)) == 4.0);

    MatrixTz az(2, 2);
    MatrixTz bz(2, 2);
    az << mpfr_float(2), mpfr_float(0), mpfr_float(0), mpfr_float(2);
    bz << mpfr_float(2), mpfr_float(0), mpfr_float(0), mpfr_float(2);
    auto cz = az * bz;
    REQUIRE(double(cz(0, 0).real()) == 4.0);
    az *= bz;
    REQUIRE(double(az(0, 0).real()) == 4.0);
  }
}
