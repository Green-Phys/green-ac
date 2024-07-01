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

#ifndef GREEN_AC_PY_NEVANLINNA_H
#define GREEN_AC_PY_NEVANLINNA_H

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <complex>

namespace green::ac::nevanlinna {

  namespace py = pybind11;
  /**
   *
   */
  py::array_t<std::complex<double>, py::array::c_style> solve(
      const py::array_t<std::complex<double> >& im_grid, const py::array_t<std::complex<double> >& grid,
      py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast>& data, double precision = 128);

}  // namespace green::ac::nevanlinna
#endif  // GREEN_AC_PY_NEVANLINNA_H