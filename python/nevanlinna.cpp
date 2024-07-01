
#include "nevanlinna.h"

#include <green/ac/nevanlinna.h>


namespace py = pybind11;

namespace green::ac::nevanlinna {

  py::array_t<std::complex<double>, py::array::c_style> solve(
      const py::array_t<std::complex<double> >& im_grid, const py::array_t<std::complex<double> >& grid,
      py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast>& data, double precision) {
    // check that frequency grids are one-dimensional arrays
    assert(grid.ndim() == 1);
    assert(im_grid.ndim() == 1);
    assert(data.ndim() > 0);
    // generate shape for output array
    std::vector<size_t> shape(data.ndim(), 0);
    shape[0] = grid.shape(0);
    if (data.ndim() > 1) std::copy(data.shape() + 1, data.shape() + data.ndim(), shape.begin() + 1);
    size_t                size = std::accumulate(shape.begin(), shape.end(), (size_t)1, std::multiplies<size_t>());
    std::complex<double>* foo  = new std::complex<double>[size];

    py::ssize_t inner_dim = std::accumulate(shape.begin() + 1, shape.end(), (size_t)1, std::multiplies<size_t>());
    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> ldata = data.reshape({data.shape(0), inner_dim});

    // Allocate local data space for analytical continuation
    ndarray::ndarray<std::complex<double>, 1> green_data(data.shape(0));
    ndarray::ndarray<std::complex<double>, 1> iwgrid(im_grid.shape(0));
    ndarray::ndarray<std::complex<double>, 1> wgrid(grid.shape(0));
    // Output array using pre-allocated dataspace
    ndarray::ndarray<std::complex<double>, 2> data_out(foo, grid.shape(0), inner_dim);
    // Fill arrays for grid points
    std::copy(im_grid.data(), im_grid.data() + im_grid.size(), iwgrid.begin());
    std::copy(grid.data(), grid.data() + grid.size(), wgrid.begin());

    green::ac::nevanlinna::nevanlinna ac(precision);
    for (size_t ind = 0; ind < inner_dim; ++ind) {
      for (size_t iw = 0; iw < data.shape(0); ++iw) {
        green_data(iw) = ldata.at(iw, ind);
      }
      ac.solve(iwgrid, green_data);
      auto out_w = ac.evaluate(wgrid);
      for (size_t iw = 0; iw < out_w.shape()[0]; ++iw) {
        data_out(iw, ind) = out_w(iw);
      }
    }

    // Create a Python object that will free the allocated
    // memory when destroyed:
    py::capsule free_when_done(foo, [](void* f) {
      std::complex<double>* foo = reinterpret_cast<std::complex<double>*>(f);
      delete[] foo;
    });

    return py::array_t<std::complex<double>, py::array::c_style>(shape,  // shape
                                                                 foo,    // the data pointer
                                                                 free_when_done);
  }

}  // namespace green::ac
