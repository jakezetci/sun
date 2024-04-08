import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import plots


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
    cmap = matplotlib.colormaps['GnBu']
    colors = cmap(np.linspace(1,0, num=18))
    for clr, j in zip(colors,density_array):

        try:
            dts = np.loadtxt(
                f'dates_{noaa_ar}_{day}_dens{j}_test2.txt', dtype=np.datetime64)
            energys = np.loadtxt(
                f'energys_11158_density={j}, day=13_test2.txt')

        except FileNotFoundError:
            continue
        ax.plot(dts, energys/(16*np.pi), '-',
                label=f'density={j}', color=clr)
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
    fig_time, ax_time = plots.config(title='time of density',
                                     xlabel='density', ylabel='time,s')
    time_spent = np.loadtxt('times_of_density2.txt')
    ax_time.plot(density_array[1:], time_spent, 'o', ms=6)
    def func_time(x, a,b):
        return (x**a) * b + 40
    coef_time, cov = curve_fit(func_time, xdata=density_array[1:], ydata=time_spent, p0=(-3, 80**3))
    print(coef_time)
    xx = np.linspace(2, 21)
    ax_time.plot(xx, (xx)**coef_time[0] * coef_time[1], label=f'y=x^{coef_time[0]:.2}', lw=2)
    fig_rel.legend()
    fig_time.legend()
    fig.legend()
    fig1.legend()
    plt.show()
