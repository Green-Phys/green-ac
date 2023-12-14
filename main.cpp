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
#include <green/ac/common_defs.h>
#include <green/ac/except.h>
#include <green/ac/nevanlinna.h>
#include <green/grids.h>
#include <green/params/params.h>
#include <green/utils/mpi_shared.h>
#include <green/utils/mpi_utils.h>
#include <mpi.h>

#include <array>

void define_parameters(green::params::params& p) {
  p.define<std::string>("input_file", "Name of the input file");
  p.define<std::string>("output_file", "Name of the output file");
  p.define<std::string>("group", "Name of the HDF5 group in the input file, that contains imaginary time data.");
  p.define<size_t>("nk", "K-point to continue, continue all if nk=0", 0);
  p.define<int>("n_omega", "Number of real frequency points", 1000);
  p.define<double>("e_min", "Smallest frequency point", -5);
  p.define<double>("e_max", "Largest frequency point", 5);
  p.define<double>("eta", "Lorentzian broadening", 0.1);
  p.define<green::ac::AC_KIND>("kind", "Type of continuation (currently only Nevanlinna is implemented)");
}

void run_nevanlinna(const green::params::params& p) {
  green::grids::transformer_t                      tr(p);
  green::ndarray::ndarray<std::complex<double>, 5> data;
  green::ndarray::ndarray<std::complex<double>, 4> data_out;
  if (!green::utils::context.global_rank) {
    green::h5pp::archive               ar(p["input_file"], "r");
    green::ndarray::ndarray<double, 1> mesh;
    ar[std::string(p["group"]) + "/data"s] >> data;
    ar[std::string(p["group"]) + "/mesh"s] >> mesh;
    if (tr.sd().repn_fermi().nts() != data.shape()[0] ||
        !std::equal(mesh.begin(), mesh.end(), tr.sd().repn_fermi().tsample().begin(),
                    [](double a, double b) { return std::abs(a - b) < 1e-12; })) {
      throw green::ac::ac_data_error("Data grid size and grid in grid_file mismatched.");
    }
    ar.close();
  }
  size_t k_shift = p["nk"];
  size_t nk = p["nk"].as<size_t>() ? 1 : data.shape()[2];
  if (size_t(p["nk"]) > data.shape()[2]) {
    throw green::ac::ac_data_error("Selected k-point number is larger than number of stored points.");
  }
  size_t ns = data.shape()[1];
  size_t nw = p["n_omega"];
  data_out.resize(std::array<size_t, 4>{nw, data.shape()[1], nk, data.shape()[3]});
  auto                              iwgrid = tr.sd().repn_fermi().wsample() * std::complex<double>(0, 1);
  std::vector<std::complex<double>> wgrid(nw);
  for (size_t iw = 0; iw < wgrid.size(); ++iw) {
    wgrid[iw] = double(p["e_min"]) + (double(p["e_max"]) - double(p["e_min"])) * double(iw) / double(wgrid.size() - 1) +
                std::complex<double>(0.0, p["eta"]);
  }
  green::ac::nevanlinna::nevanlinna ac;
  for (size_t isk = green::utils::context.global_rank; isk < nk * ns; isk += green::utils::context.global_size) {
    size_t ik = (isk % nk) + k_shift;
    size_t is = isk / nk;
    green::ndarray::ndarray<std::complex<double>, 2> inp_t(data.shape()[0], data.shape()[3]);
    green::ndarray::ndarray<std::complex<double>, 2> inp_w(tr.sd().repn_fermi().nw(), data.shape()[3]);
    for (size_t it = 0; it < data.shape()[0]; ++it) {
      for (size_t i = 0; i < data.shape()[3]; ++i) {
        inp_t(it, i) = data(it, is, ik, i, i);
      }
    }
    tr.tau_to_omega(inp_t, inp_w);
    ac.solve(iwgrid, inp_w);
    auto out_w = ac.evaluate(wgrid);
    for (size_t iw = 0; iw < out_w.shape()[0]; ++iw) {
      for (size_t i = 0; i < data.shape()[3]; ++i) {
        data_out(iw, is, ik - k_shift, i) = out_w(iw, i);
      }
    }
  }
  green::utils::allreduce(MPI_IN_PLACE, data_out.data(), data_out.size(), MPI_C_DOUBLE_COMPLEX, MPI_SUM,
                          green::utils::context.global);
  if (!green::utils::context.global_rank) {
    green::h5pp::archive ar(p["output_file"], "w");
    ar[std::string(p["group"]) + "/data"s] << data_out;
    ar[std::string(p["group"]) + "/mesh"s] << wgrid;
    ar.close();
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  std::string           name = R"(
 █▀▀█ █▀▀█ █▀▀ █▀▀ █▀▀▄
 █ ▄▄ █▄▄▀ █▀▀ █▀▀ █  █
 █▄▄█ ▀ ▀▀ ▀▀▀ ▀▀▀ ▀  ▀

 █▀▀█ █▀▀█ █▀▀▄ ▀▀█▀▀  ▀  █▀▀▄ █  █ █▀▀█ ▀▀█▀▀  ▀  █▀▀█ █▀▀▄
 █    █  █ █  █   █   ▀█▀ █  █ █  █ █▄▄█   █   ▀█▀ █  █ █  █
 █▄▄█ ▀▀▀▀ ▀  ▀   ▀   ▀▀▀ ▀  ▀  ▀▀▀ ▀  ▀   ▀   ▀▀▀ ▀▀▀▀ ▀  ▀)";
  green::params::params p(name + "\n\nGreen Analytical Continuation\n==============================================\n");
  green::grids::define_parameters(p);
  define_parameters(p);
  if (!p.parse(argc, argv)) {
    if (!green::utils::context.global_rank) p.help();
    MPI_Finalize();
    return -1;
  } else {
    if (!green::utils::context.global_rank) p.print();
  }

  try {
    switch(green::ac::AC_KIND(p["kind"])) {
      case green::ac::Nevanlinna :
        run_nevanlinna(p);
        break;
      default:
        throw std::runtime_error("Only Nevanlinna is implemented for now");
    }
  } catch (std::exception& e) {
    if (!green::utils::context.global_rank) std::cerr << e.what() << std::endl;
    MPI_Abort(green::utils::context.global, -1);
  }

  MPI_Finalize();
}