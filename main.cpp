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
  p.define<int>("nk", "K-point to continue, continue all if nk=-1", -1);
  p.define<int>("n_omega", "Number of real frequency points", 1000);
  p.define<int>("precision", "Number bit in multiprecision representation", 128);
  p.define<double>("e_min", "Smallest frequency point", -5);
  p.define<double>("e_max", "Largest frequency point", 5);
  p.define<double>("eta", "Lorentzian broadening", 0.1);
  p.define<green::ac::AC_KIND>("kind", "Type of continuation (currently only Nevanlinna is implemented)");
}


void sqrt_from_inverse(const green::ndarray::ndarray<std::complex<double>, 2>& in, green::ndarray::ndarray<std::complex<double>, 2>&& out) {
  using Matrixcd   = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using matrix_t   = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using mmatrix_t  = Eigen::Map<matrix_t>;
  using cmmatrix_t = Eigen::Map<const matrix_t>;
  size_t nso = in.shape()[0];
  Eigen::SelfAdjointEigenSolver<Matrixcd> solver(nso);
  Eigen::FullPivLU<matrix_t>              lusolver(nso, nso);
  Matrixcd S_inv = cmmatrix_t(in.data(), nso, nso);
  solver.compute(S_inv);
  mmatrix_t(out.data(), nso, nso) =
          solver.eigenvectors() * (solver.eigenvalues().cwiseSqrt().asDiagonal()) * solver.eigenvectors().adjoint();
  mmatrix_t(out.data(), nso, nso) = lusolver.compute(mmatrix_t(out.data(), nso, nso)).inverse().eval();
}

void read_nevanlinna_data(const green::params::params& p, const green::grids::transformer_t& tr,
                          green::ndarray::ndarray<std::complex<double>, 4>& data) {
  using matrix_t = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using mmatrix_t = Eigen::Map<matrix_t>;//Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  std::vector<size_t> shape(4);
  if (!green::utils::context.global_rank) {
    green::h5pp::archive               ar(p["input_file"], "r");
    green::ndarray::ndarray<double, 1> mesh;
    shape = green::h5pp::dataset_shape(ar.current_id(), std::string(p["group"]) + "/data"s);
    if (shape.size() == 4) {
      ar[std::string(p["group"]) + "/data"s] >> data;
    } else if (shape.size() == 5) {
      green::ndarray::ndarray<std::complex<double>, 5> tmp;
      ar[std::string(p["group"]) + "/data"s] >> tmp;
      size_t dim_leading = std::accumulate(tmp.shape().begin() + 1, tmp.shape().end() - 2, 1ul, std::multiplies<size_t>());
      data.resize(shape[0], shape[1], shape[2], shape[3]);
      auto tmp1 = tmp.reshape(shape[0], dim_leading, shape[3], shape[4]);
      auto tmp2 = data.reshape(shape[0], dim_leading, shape[3]);
      green::ndarray::ndarray<std::complex<double>, 3> ovlp_inv(dim_leading, shape[3], shape[4]);
      matrix_t G_orth(shape[3], shape[4]);
      ovlp_inv << -(tmp1(0) + tmp1(shape[0]-1));
      for (size_t ld = 0; ld < dim_leading; ++ld) {
        sqrt_from_inverse(ovlp_inv(ld), ovlp_inv(ld));
        for (size_t it = 0; it < shape[0]; ++it) {
          G_orth = (mmatrix_t(ovlp_inv(ld).data(), shape[3], shape[3]) * mmatrix_t(tmp1(it, ld).data(), shape[3], shape[3]) * mmatrix_t(ovlp_inv(ld).data(), shape[3], shape[3])).eval();
          for (size_t i = 0; i < shape[3]; ++i) {
            tmp2(it, ld, i) = G_orth(i,i);//tmp1(it, ld, i, i);
          }
        }
      }
    } else {
      throw green::ac::ac_data_shape_error(
          "Input data should be either 4 (diagonal orbitals)"
          " or 5 (matrix valued function) dimensional.");
    }
    ar[std::string(p["group"]) + "/mesh"s] >> mesh;
    if (tr.sd().repn_fermi().nts() != data.shape()[0] ||
        !std::equal(mesh.begin(), mesh.end(), tr.sd().repn_fermi().tsample().begin(),
                    [](double a, double b) { return std::abs(a - b) < 1e-12; })) {
      throw green::ac::ac_data_error("Data grid size and grid in grid_file mismatched.");
    }
    ar.close();
  }
  MPI_Bcast(shape.data(), 4, MPI_UNSIGNED_LONG, 0, green::utils::context.global);
  if (green::utils::context.global_rank) {
    data.resize(shape);
  }
  MPI_Bcast(data.data(), data.size(), MPI_CXX_DOUBLE_COMPLEX, 0, green::utils::context.global);
}

void run_nevanlinna(const green::params::params& p) {
  green::grids::transformer_t                      tr(p);
  green::ndarray::ndarray<std::complex<double>, 4> data;
  green::ndarray::ndarray<std::complex<double>, 4> data_out;
  read_nevanlinna_data(p, tr, data);
  size_t k_shift = p["nk"].as<int>() == -1 ? 0 : p["nk"].as<int>();
  size_t nk      = p["nk"].as<int>() == -1 ? data.shape()[2] : 1;
  size_t nao     = data.shape()[3];
  if (p["nk"].as<int>() > static_cast<int>(data.shape()[2])) {
    throw green::ac::ac_data_error("Selected k-point number is larger than number of stored points.");
  }
  size_t ns = data.shape()[1];
  size_t nw = p["n_omega"];
  data_out.resize(std::array<size_t, 4>{nw, data.shape()[1], nk, data.shape()[3]});
  auto iwgrid = tr.sd().repn_fermi().wsample() * std::complex<double>(0, 1);
  data_out.set_zero();
  green::ndarray::ndarray<std::complex<double>, 1> wgrid(nw);
  for (size_t iw = 0; iw < wgrid.size(); ++iw) {
    wgrid(iw) = double(p["e_min"]) + (double(p["e_max"]) - double(p["e_min"])) * double(iw) / double(wgrid.size()) +
                std::complex<double>(0.0, p["eta"]);
  }
  size_t dim       = (nk * ns * nao);
  double progress  = 0;
  size_t last_step = 0;
  double step_size = dim > green::utils::context.global_size ? (50.0 / (dim / green::utils::context.global_size)) : 0;
  if (!green::utils::context.global_rank) {
    std::cout << "Performing " + std::to_string(nk * ns * nao) + " continuations." << std::endl;
    std::cout << "0%|" << std::setw(55) << std::right << "|100%" << std::endl << "  |" << std::flush;
  }
  for (size_t iski = green::utils::context.global_rank; iski < dim; iski += green::utils::context.global_size) {
    size_t                                           i  = iski % nao;
    size_t                                           ik = ((iski / nao) % nk) + k_shift;
    size_t                                           is = iski / nao / nk;
    green::ndarray::ndarray<std::complex<double>, 1> inp_t(data.shape()[0]);
    green::ndarray::ndarray<std::complex<double>, 1> inp_w(tr.sd().repn_fermi().nw());
    for (size_t it = 0; it < data.shape()[0]; ++it) {
      inp_t(it) = data(it, is, ik, i);
    }
    tr.tau_to_omega(inp_t, inp_w);
    green::ac::nevanlinna::nevanlinna ac(p["precision"].as<int>());
    ac.solve(iwgrid, inp_w);
    auto out_w = ac.evaluate(wgrid);
    for (size_t iw = 0; iw < out_w.shape()[0]; ++iw) {
      data_out(iw, is, ik - k_shift, i) = out_w(iw);
    }
    if ((progress + 1) < ((iski * 100.0) / dim)) {
      if (!green::utils::context.global_rank)
        for (int step = last_step, curr = 0; step < 50 && curr < step_size; ++step, ++curr) {
          std::cout << "=" << std::flush;
          ++last_step;
        }
      progress = ((iski * 100.0) / dim);
    }
  }
  if (!green::utils::context.global_rank)
    for (int step = last_step; step < 50; ++step) std::cout << "=" << std::flush;
  if (!green::utils::context.global_rank) std::cout << "|" << std::endl;
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
  }
  if (!green::utils::context.global_rank) p.print();

  try {
    switch (p["kind"].as<green::ac::AC_KIND>()) {
      case green::ac::Nevanlinna:
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
