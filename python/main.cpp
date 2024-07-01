#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <complex>

#include "nevanlinna.h"

#define STRINGIFY(x)       #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

namespace green::ac {

  py::array_t<std::complex<double>, py::array::c_style> solve(
      const std::string& kind,
      const py::array_t<std::complex<double> >& im_grid, const py::array_t<std::complex<double> >& grid,
      py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast>& data, double precision = 128) {
    std::string _kind = kind;
    //this will make the string into upper case
    std::transform(_kind.begin(), _kind.end(), _kind.begin(), [](unsigned char c){ return std::toupper(c); });
    if ( _kind == "NEVANLINNA") {
      return nevanlinna::solve(im_grid, grid, data, precision);
    }
    throw std::runtime_error("Only kind=`Nevanlinna` continuation is supported");
  }

  PYBIND11_MODULE(_green_ac, m) {
    m.doc() = R"pbdoc(
    Green Software Package analytical continuation suit
    -----------------------
    .. currentmodule:: green_ac
    .. autosummary::
       :toctree: _generate
       solve
    )pbdoc";

    m.def("solve", &solve, R"pbdoc(
      Solve analytical continuation problem
      Arguments
      ---------
      kind      : Continuation solver, currently only 'Nevanlinna' is supported
      im_grid   : Matsubara frequency grid
      grid      : real frequency grid
      data      : Matsubara frequency data, we assume that leading dimension is frequency,
                  then continuation will be done for all inner dimensions separately
      precision : GMP precision in bits
    )pbdoc", py::arg("kind"), py::arg("im_grid"), py::arg("grid"), py::arg("data"), py::arg("precision") = 128);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
  }

}  // namespace green::ac
