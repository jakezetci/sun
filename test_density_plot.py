import matplotlib.pyplot as plt
import numpy as np


N = 16  # количество временных точек
density = 4  # пикселей на один объёмный куб
day = 13  # начальная дата
year = 2011
month = '02'
noaa_ar = 11158



if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1)
    fig_rel, ax_rel = plt.subplots(1, 1)
    density_array = np.hstack((np.linspace(1, 10, num=10-1), np.linspace(10.125, 20.125, 9)))
    print(density_array)
    ideal = np.loadtxt(
        f'energys_11158_density={density_array[1]}, day=13_test2.txt')
    error_mean = []
    for j in density_array:

        try:
            dts = np.loadtxt(
                f'dates_{noaa_ar}_{day}_dens{j}_test2.txt', dtype=np.datetime64)
            energys = np.loadtxt(
                f'energys_11158_density={j}, day=13_test2.txt')

        except FileNotFoundError:
            continue
        ax.plot(dts, energys/(16*np.pi), '-',
                label=f'density={j}')
        ideal_copy = np.copy(ideal)
        rel_error = np.abs(1-np.divide(energys, ideal))
        ideal_copy, rel_error = (np.array(t) for t in zip(*sorted(zip(ideal_copy, rel_error))))
        ax_rel.plot(ideal_copy/(16*np.pi), rel_error, 'o',
                    label=f'density={j}')
        error_mean.append(np.mean(rel_error))
        
    fig1, ax1 = plt.subplots(1,1)
    from scipy.optimize import curve_fit
    def func(x, a,b):
        return (x-2.125)**a *b
    coef, cov = curve_fit(func, xdata=density_array[1:], ydata=error_mean, p0=(2, 0.01))


    ax1.plot(density_array[1:], error_mean, label='error of density')
    ax1.plot(density_array[1:], coef[1]*(density_array[1:]-2.125)**coef[0], label=f'y=x^{coef[0]:.2}')
    fig_rel.legend()
    fig.legend()
    fig1.legend()
    plt.show()
