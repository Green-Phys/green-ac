import importlib.util as ilu
import numpy as np
import os

import concurrent.futures
from multiprocessing import cpu_count, Process

has_mpi = False

if ilu.find_spec("mpi4py") is not None:
    from mpi4py import MPI
    has_mpi = True

from ._green_ac import __version__
from ._green_ac import solve as _solve


def solve_mpi(comm, kind, mgrid, grid, data, prec=128):
    '''
    Solve analytical continuation problem using MPI parallelization. MPI4Py should be available and properly built
    '''
    if not has_mpi:
        raise RuntimeError("MPI is not found")
    lin_s = data.size//data.shape[0]
    data_rs = data.reshape([data.shape[0], lin_s])
    out_data_rs = np.zeros([grid.shape[0], lin_s], dtype=np.complex128)
    rank = comm.Get_rank()
    size = comm.Get_size()
    for i in range(rank, lin_s, size):
        out_data_rs[:,i] = _solve(kind, mgrid, grid, data_rs[:,i], prec)
    comm.Allreduce(MPI.IN_PLACE, out_data_rs.view(np.float64), op=MPI.SUM)
    out_data = out_data_rs.reshape(grid.shape + data.shape[1:])
    return out_data

def _solve_conc(kind, mgrid, grid, data_rs, prec):
    return _solve(kind, mgrid, grid, data_rs, prec)

def solve(kind, mgrid, grid, data, prec=128):
    '''
    Solve analytical continuation using Python built-in concurrent parallelism
    '''
    from multiprocessing import cpu_count, Process
    # Find the number of cpus to use for parallelization of analtyic continuation
    slurm_env_var = 'SLURM_JOB_CPUS_PER_NODE'
    phys_ncpu = cpu_count()
    _ncpus = int(os.environ[slurm_env_var]) if slurm_env_var in os.environ else int(phys_ncpu // 2)
    # linearize data for easy chunking
    lin_s = data.size//data.shape[0]
    data_rs = data.reshape([data.shape[0], lin_s])
    out_data_rs = np.zeros([grid.shape[0], lin_s], dtype=np.complex128)
    fs = []
    chunks = np.array_split(range(lin_s), _ncpus)
    # create a process pool
    with concurrent.futures.ProcessPoolExecutor(_ncpus) as executor:
        for i in chunks:
            if len(i) == 0:
                continue
            print(i)
            f = executor.submit(_solve_conc, kind, mgrid, grid, data_rs[:, i], prec)
            fs.append(f)
    # wait for completetion
    concurrent.futures.wait(fs)
    # get data from each process
    i = 0
    for c in (chunks):
        if len(c) == 0:
            continue
        out_data_rs[:,c] = fs[i].result()
        i += 1
    out_data = out_data_rs.reshape(grid.shape + data.shape[1:])
    return out_data
