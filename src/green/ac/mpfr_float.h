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

#ifndef AC_MPFR_FLOAT_H
#define AC_MPFR_FLOAT_H

#include <Eigen/Core>
#include <complex>
#include <type_traits>

#include "mpfr.h"
#include "mpfr_float.h"

namespace green::ac {

  // forward declaration
  struct mpfr_float;

  namespace type_traits {
    template <typename T>
    using convertible_t = std::enable_if_t<!std::is_same_v<std::decay_t<T>, mpfr_float> && std::is_arithmetic_v<T>, mpfr_float>;
    template <typename T>
    constexpr bool is_scalar = std::is_same_v<std::decay_t<T>, mpfr_float> || std::is_arithmetic_v<T>;
    template <typename T>
    using scalar_t = std::enable_if_t<is_scalar<T>, T>;

  }  // namespace type_traits

  inline static int precision = 100;

  struct mpfr_float {
    mpfr_float() : mpfr_float(0.0) {}
    mpfr_float(double rhs) {
      mpfr_init2(val, precision);
      mpfr_set_d(val, rhs, MPFR_RNDD);
    }

    ~mpfr_float() { mpfr_clear(val); }
    /*
     * Copy/Move constructors
     */
    mpfr_float(const mpfr_float& rhs) {
      mpfr_init2(val, precision);
      mpfr_set(val, rhs.val, MPFR_RNDD);
    }
    mpfr_float(mpfr_float&& rhs) noexcept {
      val->_mpfr_d = nullptr;
      mpfr_swap(val, rhs.val);
    }

    /*
     * Copy/Move assignment
     */
    inline mpfr_float& operator=(const mpfr_float& rhs) {
      if (&rhs == this) return *this;
      mpfr_set(val, rhs.val, MPFR_RNDD);
      return *this;
    }
    inline mpfr_float& operator=(mpfr_float&& rhs) noexcept {
      if (&rhs == this) return *this;
      mpfr_swap(val, rhs.val);
      return *this;
    }

    inline            operator double() const { return mpfr_get_d(val, MPFR_RNDD); }

    mpfr_t            val{};

    inline mpfr_float operator-() const {
      mpfr_float res(*this);
      return res *= -1;
    }

#define MPFR_MATH_OP(OP, MPFR_OP)                      \
  template <typename T>                                \
  inline mpfr_float& operator OP(T && rhs) {           \
    MPFR_OP(val, val, mpfr_float(rhs).val, MPFR_RNDD); \
    return *this;                                      \
  }

    MPFR_MATH_OP(+=, mpfr_add)
    MPFR_MATH_OP(-=, mpfr_sub)
    MPFR_MATH_OP(*=, mpfr_mul)
    MPFR_MATH_OP(/=, mpfr_div)

    template <typename T>
    T convert_to() const {
      return T(mpfr_get_d(val, MPFR_RNDD));
    }
  };

#define MPFR_BIN_MATH_OP(OP, MPFR_OP)                                                  \
  template <typename T>                                                                \
  inline type_traits::convertible_t<T> operator OP(const mpfr_float & lhs, T && rhs) { \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, lhs.val, mpfr_float(rhs).val, MPFR_RNDD);                         \
    return res;                                                                        \
  }                                                                                    \
  template <typename T>                                                                \
  inline type_traits::convertible_t<T> operator OP(T && lhs, const mpfr_float & rhs) { \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, mpfr_float(lhs).val, rhs.val, MPFR_RNDD);                         \
    return res;                                                                        \
  }                                                                                    \
  template <typename T>                                                                \
  inline type_traits::convertible_t<T> operator OP(mpfr_float && lhs, T && rhs) {      \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, lhs.val, mpfr_float(rhs).val, MPFR_RNDD);                         \
    return res;                                                                        \
  }                                                                                    \
  template <typename T>                                                                \
  inline type_traits::convertible_t<T> operator OP(T && lhs, mpfr_float && rhs) {      \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, mpfr_float(lhs).val, rhs.val, MPFR_RNDD);                         \
    return res;                                                                        \
  }                                                                                    \
  inline mpfr_float operator OP(const mpfr_float& lhs, const mpfr_float& rhs) {        \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, lhs.val, mpfr_float(rhs).val, MPFR_RNDD);                         \
    return res;                                                                        \
  }                                                                                    \
  inline mpfr_float operator OP(mpfr_float&& lhs, mpfr_float&& rhs) {                  \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, mpfr_float(lhs).val, rhs.val, MPFR_RNDD);                         \
    return res;                                                                        \
  }                                                                                    \
  inline mpfr_float operator OP(const mpfr_float& lhs, mpfr_float&& rhs) {             \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, lhs.val, mpfr_float(rhs).val, MPFR_RNDD);                         \
    return res;                                                                        \
  }                                                                                    \
  inline mpfr_float operator OP(mpfr_float&& lhs, const mpfr_float& rhs) {             \
    mpfr_float res(0.0);                                                               \
    MPFR_OP(res.val, mpfr_float(lhs).val, rhs.val, MPFR_RNDD);                         \
    return res;                                                                        \
  }

  MPFR_BIN_MATH_OP(+, mpfr_add)
  MPFR_BIN_MATH_OP(-, mpfr_sub)
  MPFR_BIN_MATH_OP(*, mpfr_mul)
  MPFR_BIN_MATH_OP(/, mpfr_div)

#define MPFR_MATH_FUN_OP(FUN, MPFR_FUN)          \
  inline mpfr_float FUN(const mpfr_float& rhs) { \
    mpfr_float res(0.0);                         \
    MPFR_FUN(res.val, rhs.val, MPFR_RNDD);       \
    return res;                                  \
  }

  MPFR_MATH_FUN_OP(sqrt, mpfr_sqrt)
  MPFR_MATH_FUN_OP(cos, mpfr_cos)
  MPFR_MATH_FUN_OP(sin, mpfr_sin)

  inline mpfr_float atan2(const mpfr_float& x, const mpfr_float& y) {
    mpfr_float res(0.0);
    mpfr_atan2(res.val, x.val, y.val, MPFR_RNDD);
    return res;
  }

}  // namespace green::ac

namespace std {
  // std::complex explicit template specialization for mpfr_float to overcome libc++ limitations.
  template <>
  class complex<green::ac::mpfr_float> {
    using real_t = green::ac::mpfr_float;

  public:
    complex() : _real(real_t(0)), _imag(real_t(0)) {}
    complex(const real_t& real, const real_t& imag) : _real(real), _imag(imag) {}
    complex(const real_t& real) : _real(real), _imag(real_t(0)) {}

    template <typename S>
    complex(const std::enable_if_t<std::is_arithmetic_v<S>, S>& real) : _real(real_t(real)), _imag(real_t(0)) {}

    template <typename S>
    complex(const std::complex<S>& cplx) :
        _real(green::ac::type_traits::convertible_t<S>(cplx.real())), _imag(real_t(cplx.imag())) {}

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator=(const S& rhs) {
      _real = real_t(rhs);
      _imag = real_t(0);
      return *this;
    }
    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator=(const std::complex<S>& rhs) {
      _real = real_t(rhs.real());
      _imag = real_t(rhs.imag());
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator-=(const complex<S>& rhs) {
      _real -= real_t(rhs.real());
      _imag -= real_t(rhs.imag());
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator+=(const complex<S>& rhs) {
      _real += real_t(rhs.real());
      _imag += real_t(rhs.imag());
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator*=(const complex<S>& rhs) {
      *this = *this * rhs;
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator/=(const complex<S>& rhs) {
      *this = *this / rhs;
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator-=(const S& rhs) {
      _real -= real_t(rhs);
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator+=(const S& rhs) {
      _real += real_t(rhs);
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator*=(const S& rhs) {
      _real *= real_t(rhs);
      _imag *= real_t(rhs);
      return *this;
    }

    template <typename S>
    std::enable_if_t<green::ac::type_traits::is_scalar<S>, complex>& operator/=(const S& rhs) {
      _real /= real_t(rhs);
      _imag /= real_t(rhs);
      return *this;
    }

    complex operator-() const { return complex{-_real, -_imag}; }

  private:
    real_t _real;
    real_t _imag;

  public:
    [[nodiscard]] real_t real() const { return _real; }
    [[nodiscard]] real_t imag() const { return _imag; }

    void                 real(const real_t& re) { _real = re; }
    void                 imag(const real_t& im) { _imag = im; }
  };

  template <>
  inline bool operator==(const complex<green::ac::mpfr_float>& x, const complex<green::ac::mpfr_float>& y) {
    return x.real() == y.real() && x.imag() == y.imag();
  }

  template <typename S>
  inline bool operator==(const complex<green::ac::mpfr_float>& x, const green::ac::type_traits::convertible_t<S>& y) {
    return x.real() == y && x.imag() == 0;
  }

  template <typename S>
  inline bool operator==(const green::ac::type_traits::convertible_t<S>& y, const complex<green::ac::mpfr_float>& x) {
    return x.real() == y && x.imag() == 0;
  }

  inline auto operator*(const complex<green::ac::mpfr_float>& x, const complex<green::ac::mpfr_float>& y) {
    green::ac::mpfr_float a = x.real();
    green::ac::mpfr_float b = x.imag();
    green::ac::mpfr_float c = y.real();
    green::ac::mpfr_float d = y.imag();
    return std::complex<green::ac::mpfr_float>((a * c - b * d), (a * d + b * c));
  }

  inline auto operator/(const complex<green::ac::mpfr_float>& x, const complex<green::ac::mpfr_float>& y) {
    green::ac::mpfr_float a     = x.real();
    green::ac::mpfr_float b     = x.imag();
    green::ac::mpfr_float c     = y.real();
    green::ac::mpfr_float d     = y.imag();
    green::ac::mpfr_float denom = c * c + d * d;
    green::ac::mpfr_float r     = (a * c + b * d) / denom;
    green::ac::mpfr_float i     = (b * c - a * d) / denom;
    return typename std::complex<green::ac::mpfr_float>(r, i);
  }

  inline auto operator+(const complex<green::ac::mpfr_float>& x, const complex<green::ac::mpfr_float>& y) {
    green::ac::mpfr_float a = x.real();
    green::ac::mpfr_float b = x.imag();
    green::ac::mpfr_float c = y.real();
    green::ac::mpfr_float d = y.imag();
    return typename std::complex<green::ac::mpfr_float>((a + c), (b + d));
  }

  inline auto operator-(const complex<green::ac::mpfr_float>& x, const complex<green::ac::mpfr_float>& y) {
    green::ac::mpfr_float a = x.real();
    green::ac::mpfr_float b = x.imag();
    green::ac::mpfr_float c = y.real();
    green::ac::mpfr_float d = y.imag();
    return std::complex<green::ac::mpfr_float>((a - c), (b - d));
  }

  inline green::ac::mpfr_float abs(const complex<green::ac::mpfr_float>& x) {
    return green::ac::sqrt(x.real() * x.real() + x.imag() * x.imag());
  }

  inline complex<green::ac::mpfr_float> conj(const complex<green::ac::mpfr_float>& x) {
    return complex<green::ac::mpfr_float>(x.real(), -x.imag());
  }

  inline complex<green::ac::mpfr_float> polar(const green::ac::mpfr_float& rho, const green::ac::mpfr_float& arg) {
    green::ac::mpfr_float x = green::ac::cos(arg);
    green::ac::mpfr_float y = green::ac::sin(arg);
    x *= rho;
    y *= rho;
    return complex<green::ac::mpfr_float>(x, y);
  }

  inline green::ac::mpfr_float arg(const complex<green::ac::mpfr_float>& x) {
    green::ac::mpfr_float arg_;
    arg_ = green::ac::atan2(x.imag(), x.real());
    return arg_;
  }

  inline complex<green::ac::mpfr_float> sqrt(const complex<green::ac::mpfr_float>& x) {
    return std::polar(green::ac::sqrt(std::abs(x)), std::arg(x) / green::ac::mpfr_float(2));
  }

  inline green::ac::mpfr_float real(const complex<green::ac::mpfr_float>& x) { return x.real(); }

  inline green::ac::mpfr_float imag(const complex<green::ac::mpfr_float>& x) { return x.imag(); }
}  // namespace std

namespace Eigen {

  template <>
  struct NumTraits<std::complex<green::ac::mpfr_float>> :
      NumTraits<std::complex<double>>  // permits to get the epsilon, dummy_precision, lowest, highest functions
  {
    typedef green::ac::mpfr_float Real;
    typedef green::ac::mpfr_float NonInteger;
    typedef green::ac::mpfr_float Nested;

    enum { IsComplex = 1, IsInteger = 0, IsSigned = 1, RequireInitialization = 1, ReadCost = 1, AddCost = 3, MulCost = 3 };
  };

}  // namespace Eigen

#endif  // AC_MPFR_FLOAT_H
