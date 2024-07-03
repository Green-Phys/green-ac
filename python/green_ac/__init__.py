import concurrent.futures
import functools
import importlib.util as ilu
import multiprocessing
import os

import numpy as np


def default_kwargs(**defaultKwargs):
    def actual_decorator(fn):
        @functools.wraps(fn)
        def g(*args, **kwargs):
            defaultKwargs.update(kwargs)
            return fn(*args, **defaultKwargs)

        return g

    return actual_decorator


has_mpi = False

if ilu.find_spec("mpi4py") is not None:
    from mpi4py import MPI

    has_mpi = True

from ._green_ac import solve as _solve


def _solve_nevanlinna_conc(kind, mgrid, grid, data_rs, prec):
    """
    This is a wrapper around C++ function to make it work with ProcessPoolExecutor
    """
    return _solve(kind, mgrid, grid, data_rs, prec)


def _solve_nevanlinna_mpi():
    pass


@default_kwargs(eta=0.01, prec=128)
def _solve_nevanlinna(mgrid, grid, data, **kwargs):
    eta = kwargs["eta"]
    prec = kwargs["prec"]
    from multiprocessing import cpu_count
    # Find the number of cpus to use for parallelization of analytic continuation
    slurm_env_var = 'SLURM_JOB_CPUS_PER_NODE'
    phys_ncpu = cpu_count()
    _ncpus = int(os.environ[slurm_env_var]) if slurm_env_var in os.environ else int(phys_ncpu // 2)
    # linearize data for easy chunking
    lin_s = data.size // data.shape[0]
    data_rs = data.reshape([data.shape[0], lin_s])
    out_data_rs = np.zeros([grid.shape[0], lin_s], dtype=np.complex128)
    fs = []
    chunks = np.array_split(range(lin_s), _ncpus)
    # create a process pool
    with concurrent.futures.ProcessPoolExecutor(_ncpus, mp_context=multiprocessing.get_context('forkserver')) as executor:
        for i in chunks:
            if len(i) == 0:
                continue
            f = executor.submit(_solve_nevanlinna_conc, "Nevanlinna", mgrid * 1.j, grid + eta * 1.j, data_rs[:, i],
                                prec)
            fs.append(f)
    # wait for completion
    concurrent.futures.wait(fs)
    # get data from each process
    i = 0
    for c in (chunks):
        if len(c) == 0:
            continue
        out_data_rs[:, c] = fs[i].result()
        i += 1
    out_data = out_data_rs.reshape(grid.shape + data.shape[1:])
    return out_data


solvers = {"NEVANLINNA": _solve_nevanlinna}
solvers_mpi = {"NEVANLINNA": _solve_nevanlinna_mpi}


@default_kwargs(eta=0.01, prec=128)
def solve_mpi(comm, kind, mgrid, grid, data, **kwargs):
    """
    Solve analytical continuation problem using MPI parallelization
    MPI4Py should be available and properly built

    :param comm: MPI communicator
    :param kind: type of continuation (should be string)
    :param mgrid: Matsubara frequency points (real one dimensional numpy array)
    :param grid: real frequency points (real one dimensional numpy array)
    :param data: Matsubara frequency Green's function (complex N-dimensional numpy array)
    :param kwargs: additional solver specific parameters such as broadening `eta` or precision
    :return: complex N-dimensional numpy array with real frequency Green's function
    """
    if not has_mpi:
        raise RuntimeError("MPI is not found")
    prec = kwargs["prec"]
    eta = kwargs["eta"]
    lin_s = data.size // data.shape[0]
    data_rs = data.reshape([data.shape[0], lin_s])
    out_data_rs = np.zeros([grid.shape[0], lin_s], dtype=np.complex128)
    rank = comm.Get_rank()
    size = comm.Get_size()
    for i in range(rank, lin_s, size):
        out_data_rs[:, i] = _solve(kind, mgrid * 1.j, grid + eta * 1.j, data_rs[:, i], prec)
    comm.Allreduce(MPI.IN_PLACE, out_data_rs.view(np.float64), op=MPI.SUM)
    out_data = out_data_rs.reshape(grid.shape + data.shape[1:])
    return out_data


def solve(kind, mgrid, grid, data, **kwargs):
    """Solve analytical continuation using Python built-in concurrent parallelism

    :param kind: type of continuation (should be string)
    :param mgrid: Matsubara frequency points (real one dimensional numpy array)
    :param grid: real frequency points (real one dimensional numpy array)
    :param data: Matsubara frequency Green's function (complex N-dimensional numpy array)
    :param kwargs: additional solver specific parameters such as broadening `eta` or precision
    :return: complex N-dimensional numpy array with real frequency Green's function
    """
    if not kind.upper() in solvers.keys():
        raise RuntimeError(f"Unknown continuation kind {kind}, currently only {solvers.keys()} are supported.")
    return solvers[kind.upper()](mgrid, grid, data, **kwargs)
