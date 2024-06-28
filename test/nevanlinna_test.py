import pytest

import green_ac
import numpy as np


def test_two_peaks():
    imgrid = 1.j * (2*np.linspace(0,100,101) + 1) * np.pi/ 2
    grid = np.linspace(-2,2,1001) + 0.01j
    data = 0.5*(1/(imgrid + 0.5) + 1/(imgrid - 0.5))
    data_test = 0.5*(1/(grid + 0.5) + 1/(grid - 0.5))

    with pytest.raises(Exception) as e:
        green_ac.solve("Nevanlina", imgrid, grid, data)
    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data)

    assert np.allclose(data_test, data_out)

def test_two_orbitals():
    imgrid = 1.j * (2*np.linspace(0,100,101) + 1) * np.pi/ 2
    grid = np.linspace(-2,2,1001) + 0.01j
    data = np.zeros(imgrid.shape +(2,), dtype=np.complex128)
    data[:, 0] = 0.5*(1/(imgrid + 0.5) + 1/(imgrid - 0.5))
    data[:, 1] = 0.5*(1/(imgrid + 0.2) + 1/(imgrid - 0.2))
    data_test_1 = 0.5*(1/(grid + 0.5) + 1/(grid - 0.5))
    data_test_2 = 0.5*(1/(grid + 0.2) + 1/(grid - 0.2))

    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data)

    np.savetxt("Aw.txt", np.column_stack([grid.real, data_out[:, 0].imag, data_out[:,1].imag]))
    assert data_out.shape == (grid.shape[0], 2)

    assert np.allclose(data_test_1, data_out[:,0])
    assert np.allclose(data_test_2, data_out[:,1])


def test_two_orbitals_momenta():
    imgrid = 1.j * (2*np.linspace(0,100,101) + 1) * np.pi/ 2
    grid = np.linspace(-2,2,1001) + 0.01j
    data = np.zeros(imgrid.shape +(10,2,), dtype=np.complex128)
    for k in range(10):
      data[:, k, 0] = 0.5*(1/(imgrid + 0.5) + 1/(imgrid - 0.5))
      data[:, k, 1] = 0.5*(1/(imgrid + 0.2) + 1/(imgrid - 0.2))
    data_test_1 = 0.5*(1/(grid + 0.5) + 1/(grid - 0.5))
    data_test_2 = 0.5*(1/(grid + 0.2) + 1/(grid - 0.2))

    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data)

    np.savetxt("Aw.txt", np.column_stack([grid.real, data_out[:, 0].imag, data_out[:,1].imag]))
    assert data_out.shape == (grid.shape[0], 10, 2)

    assert np.allclose(data_test_1, data_out[:,0,0])
    assert np.allclose(data_test_2, data_out[:,1,1])

