import numpy as np
import lib
import computing
import cpp_module as cpp
import time
import matplotlib.pyplot as plt
import plots

if __name__ == "__main__":
    latlim = (10, 20)
    lonlim = (10, 20)

    r_lim = (696300*1e3, 1.0436*696300*1e3)
    N_2d = 4
    grid = lib.create_grid(latlim, lonlim, 4)
    m = np.array([0, 1, 0])
    grid = computing.model_grid(
        grid, m, np.zeros(3), vector=False, returnobj=True)

    N_array = np.geomspace(10, 10000, num=50)
    energys = np.zeros_like(N_array)
    values, points, areas = grid.values, grid.values, grid.area

    tic = time.perf_counter()
    timestamp = 10
    for i, N in enumerate(N_array):
        toc = time.perf_counter()
        gridvol, vol = lib.create_3Dgrid_sph(latlim, lonlim, r_lim, int(N))
        energy = 0
        for xyz in grid.xyz:
            B = cpp.b_comp(xyz, values, points, areas)
            energy = energy + (np.linalg.norm(B)**2 * vol)
        energys[i] = energy / (8*np.pi)
        if i % timestamp == 0:
            toc = time.perf_counter()
            print(
                f"values {i-timestamp} - {i} done in {toc - tic:0.2f} seconds")
            tic = time.perf_counter()

    np.savetxt('energy_of_density.txt', energys)

    fig, ax = plots.config(xlabel='Density', ylabel='Energy', logscalex=True)
    ax.plot(N_array, energys, 'o')
