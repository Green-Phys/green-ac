[![GitHub license](https://img.shields.io/github/license/Green-Phys/green-ac?cacheSeconds=3600&color=informational&label=License)](./LICENSE)
[![GitHub license](https://img.shields.io/badge/C%2B%2B-17-blue)](https://en.cppreference.com/w/cpp/compiler_support/17)

![grids](https://github.com/Green-Phys/green-ac/actions/workflows/test.yaml/badge.svg)
[![codecov](https://codecov.io/github/Green-Phys/green-ac/graph/badge.svg?token=MWM2Y2HR8X)](https://codecov.io/github/Green-Phys/green-ac)


```
 █▀▀█ █▀▀█ █▀▀ █▀▀ █▀▀▄
 █ ▄▄ █▄▄▀ █▀▀ █▀▀ █  █
 █▄▄█ ▀ ▀▀ ▀▀▀ ▀▀▀ ▀  ▀

 █▀▀█ █▀▀█ █▀▀▄ ▀▀█▀▀  ▀  █▀▀▄ █  █ █▀▀█ ▀▀█▀▀  ▀  █▀▀█ █▀▀▄
 █    █  █ █  █   █   ▀█▀ █  █ █  █ █▄▄█   █   ▀█▀ █  █ █  █
 █▄▄█ ▀▀▀▀ ▀  ▀   ▀   ▀▀▀ ▀  ▀  ▀▀▀ ▀  ▀   ▀   ▀▀▀ ▀▀▀▀ ▀  ▀
```
***

`Green/Continuation` is an Analytical continuation toolkit for Green Software Package

# Installation

`Green/Continuation` comes in two forms

  - C++ application
  - Python package

## Dependencies

`Green/Continuation` has the following required external dependencies
  - HDF5 library version >= 1.10.2
  - Message Passing Interface >= 3.1 (for C++ application)
  - Eigen3 library >= 3.4.0
  - GNU Multiprecision library
  - pybind11 (optional to build python wrapper)

To build `Green/Continuation` CMake version 3.18 or above is required

## Build and Install C++ application

The following example will build, test and install `Green/Continuation` to `/path/to/weakcoupling/install/dir` directory.

```ShellSession
$ git clone https://github.com/Green-Phys/green-ac
$ cd green-ac
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..
$ make
$ make test
$ make install
```

## Installation of Python package

We provide pre-built binaries for major Linux distributions and recent MacOS version via `pip`.
For installation using `pip `simply type `pip install green-ac`.
If pre-built binaries can not be used, package will be built from sources.

## Usage

### C++ application

After the `Green/Continuation` is built and installed, spectral function could be obtained by calling

```ShellSession
<install dir>/bin/ac.exe  --BETA <Inverse temperature> --grid_file <grid file>    \
   --input_file <input file> --output_file <output file> --group <HDF5 group with data> \
   --e_min -5.0 --e_max 5.0 --n_omega 4000 --eta 0.01 \
   --kind Nevanlinna
```

  - `BETA`  -- the inverse temperature that was used to obtain results
  - `grid_file` -- name of the grid file that was used to obtain results
  - `input_file` -- name of the file that contains imaginary time data
  - `output_file` -- name of the file to store results of the continuation
  - `group` -- name of the group that contains `data` and `mesh` datasets with imaginary time data and grid
  - `e_min` -- lowest frequency on the real axis
  - `e_max` -- largest frequency on the real axis
  - `n_omega` -- number of frequency points on the real axis
  - `eta` -- broadening parameter
  - `kind` -- type of continuation to be used (current version only supports `Nevanlinna`)

After the completetion, results will be stored in `group` HDF5 group in the `output_file` file.


### Python

In addition to C++ application, `Green/Continuation` provides a convinient Python package that can work directly with
`numpy` arrays. It supports two types of parallelism, using `ProcessPoolExecutor` from `concurrent.futures` and
MPI parallelization using `mpi4py` library.

To use `Green/Continuation` simply import it in your script as

```Python
import green_ac
```

and call `solve` function with the following parameters:

  - Type of the continuation (currently we only provide `Nevanlinna`)
  - Matsubara frequency grid
  - Real frequency grid
  - Data in Matsubara frequency domain
  - Precision

Here is an example how to obtain real frequency Green's function for a simple two-pole non-interacting Green's function with precision at least 512 bits:


```Python
imgrid = (2*np.linspace(-50,49,100) + 1) * 1.j *np.pi/ 2
grid = np.linspace(-2,2,1001) + 0.01j
data = 0.5*(1/(imgrid + 0.5) + 1/(imgrid - 0.5))
data_out = green_ac.solve("Nevanlinna", imgrid, grid, data, 512)
```

Here we have `100` positive and negative Matsubara frequencies, define real frequency grid to have `1000` points and to be from `-2` to `2` with broadening parameter `0.01`,
We define Green's function on Matsubara grid with poles at `-0.5` and `0.5`. Real frequency data will be stored in `data_out` array.


# Acknowledgements

This work is supported by National Science Foundation under the award OCA-2310582
