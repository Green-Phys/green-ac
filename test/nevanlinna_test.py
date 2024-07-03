import os
import pytest

import green_ac
import numpy as np


def test_two_peaks():
    imgrid = (2 * np.linspace(-50, 49, 100) + 1) * np.pi / 2
    grid = np.linspace(-2, 2, 1001)
    eta = 0.01
    data = 0.5 * (1 / (1.j * imgrid + 0.5) + 1 / (1.j * imgrid - 0.5))
    data_test = 0.5 * (1 / (grid + eta * 1.j + 0.5) + 1 / (grid + eta * 1.j - 0.5))

    with pytest.raises(Exception) as e:
        green_ac.solve("Nevanlina", imgrid, grid, data, eta=eta)
    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data, eta=eta)

    assert np.allclose(data_test, data_out)


def test_default_eta():
    imgrid = (2 * np.linspace(-50, 49, 100) + 1) * np.pi / 2
    grid = np.linspace(-2, 2, 1001)
    eta = 0.01
    data = 0.5 * (1 / (1.j * imgrid + 0.5) + 1 / (1.j * imgrid - 0.5))
    data_test = 0.5 * (1 / (grid + eta * 1.j + 0.5) + 1 / (grid + eta * 1.j - 0.5))

    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data)

    assert np.allclose(data_test, data_out)


def test_two_orbitals():
    imgrid = (2 * np.linspace(0, 100, 101) + 1) * np.pi / 2
    grid = np.linspace(-2, 2, 1001)
    eta = 0.01
    data = np.zeros(imgrid.shape + (2,), dtype=np.complex128)
    data[:, 0] = 0.5 * (1 / (1.j * imgrid + 0.5) + 1 / (1.j * imgrid - 0.5))
    data[:, 1] = 0.5 * (1 / (1.j * imgrid + 0.2) + 1 / (1.j * imgrid - 0.2))
    data_test_1 = 0.5 * (1 / (grid + eta * 1.j + 0.5) + 1 / (grid + eta * 1.j - 0.5))
    data_test_2 = 0.5 * (1 / (grid + eta * 1.j + 0.2) + 1 / (grid + eta * 1.j - 0.2))

    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data, eta=eta)

    assert data_out.shape == (grid.shape[0], 2)

    assert np.allclose(data_test_1, data_out[:, 0])
    assert np.allclose(data_test_2, data_out[:, 1])


def test_two_orbitals_momenta():
    imgrid = (2 * np.linspace(0, 100, 101) + 1) * np.pi / 2
    grid = np.linspace(-2, 2, 1001)
    eta = 0.01
    data = np.zeros(imgrid.shape + (10, 2,), dtype=np.complex128)
    for k in range(10):
        data[:, k, 0] = 0.5 * (1 / (1.j * imgrid + 0.5) + 1 / (1.j * imgrid - 0.5))
        data[:, k, 1] = 0.5 * (1 / (1.j * imgrid + 0.2) + 1 / (1.j * imgrid - 0.2))
    data_test_1 = 0.5 * (1 / (grid + eta * 1.j + 0.5) + 1 / (grid + eta * 1.j - 0.5))
    data_test_2 = 0.5 * (1 / (grid + eta * 1.j + 0.2) + 1 / (grid + eta * 1.j - 0.2))

    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data, eta=eta)

    assert data_out.shape == (grid.shape[0], 10, 2)

    assert np.allclose(data_test_1, data_out[:, 0, 0])
    assert np.allclose(data_test_2, data_out[:, 1, 1])


def test_slicing():
    imgrid = (2 * np.linspace(0, 100, 101) + 1) * np.pi / 2
    grid = np.linspace(-2, 2, 1001)
    eta = 0.01
    N = 10
    data = np.zeros(imgrid.shape + (N, 2,), dtype=np.complex128)
    data_test = np.zeros(grid.shape + (N, 2,), dtype=np.complex128)
    eps_k = np.cos(np.linspace(0, 2 * np.pi, N, endpoint=False))
    for k in range(10):
        data[:, k, 0] = 0.5 * (1 / (1.j * imgrid + 0.5 - eps_k[k]) + 1 / (1.j * imgrid - 0.5 - eps_k[k]))
        data[:, k, 1] = 0.5 * (1 / (1.j * imgrid + 0.2 - eps_k[k]) + 1 / (1.j * imgrid - 0.2 - eps_k[k]))
        data_test[:, k, 0] = 0.5 * (1 / (grid + eta * 1.j + 0.5 - eps_k[k]) + 1 / (grid + eta * 1.j - 0.5 - eps_k[k]))
        data_test[:, k, 1] = 0.5 * (1 / (grid + eta * 1.j + 0.2 - eps_k[k]) + 1 / (grid + eta * 1.j - 0.2 - eps_k[k]))
    data2 = np.copy(data)

    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data[:, 2, :], eta=eta)
    assert data_out.shape == (grid.shape[0], 2)

    assert np.allclose(data_test[:, 2, :], data_out[:, :])
    assert np.allclose(data, data2)

    data_out = green_ac.solve("Nevanlinna", imgrid, grid, data[:, 2:4, :], eta=eta)

    assert data_out.shape == (grid.shape[0], 2, 2)

    assert np.allclose(data_test[:, 3, 0], data_out[:, 1, 0])
    assert np.allclose(data_test[:, 2, 1], data_out[:, 0, 1])
    assert np.allclose(data, data2)


def test_from_file():
    data_dir = os.path.abspath(os.path.dirname(__file__)) + '/data'
    Xiw = np.loadtxt(data_dir + "/X_iw.txt")
    Xw = np.loadtxt(data_dir + "/X_w.txt")
    mgrid = Xiw[:, 0]
    grid = Xw[:, 0]
    data = Xiw[:, 1] + 1.j * Xiw[:, 2]
    test_data = Xw[:, 1]
    data_out = green_ac.solve("Nevanlinna", mgrid, grid, data, eta=0.01)
    np.allclose(test_data, -data_out.imag / np.pi, atol=1e-3)
