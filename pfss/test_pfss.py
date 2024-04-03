import pfsspy
import pfsspy.utils

import os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

if __name__ == "__main__":
    time = '2011/02/13'

    # time = a.Time('2010/01/01', '2010/01/02')
    series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')
    crot = a.jsoc.PrimeKey('CAR_ROT',
                           int(sunpy.coordinates.sun.carrington_rotation_number(time)))
    print(sunpy.coordinates.sun.carrington_rotation_number(time))
    result = Fido.search(series, crot,
                         a.jsoc.Notify('rrzhdanov@edu.hse.ru'))
    files = Fido.fetch(result)
    hmi_map = sunpy.map.Map(files[0])
    print('Data shape: ', hmi_map.data.shape)
    fig = plt.figure(figsize=(12, 5))
    ax = plt.subplot(projection=hmi_map)
    im = hmi_map.plot(axes=ax)

    ax.coords[0].set_axislabel("Carrington Longitude [deg]")
    ax.coords[1].set_axislabel("Latitude [deg]")

    ax.coords.grid(color='black', alpha=0.6, linestyle='dotted', linewidth=0.5)

    cb = plt.colorbar(im, fraction=0.019, pad=0.1)
    cb.set_label(f"Radial magnetic field [{hmi_map.unit}]")

    # In order to make the x-axis ticks show, the bottom y-limit has to be adjusted slightly
    ax.set_ylim(bottom=0)
    ax.set_title(f"{hmi_map.meta['content']},\n"
                 f"Carrington rotation {hmi_map.meta['CAR_ROT']}")

    plt.show()
    hmi_map = hmi_map.resample([360, 180] * u.pix)
    nrho = 35
    rss = 1.08
    pfss_in = pfsspy.Input(hmi_map, nrho, rss)
    pfss_out = pfsspy.pfss(pfss_in)
    ss_br = pfss_out.source_surface_br
    # Create the figure and axes
    fig = plt.figure()
    ax = plt.subplot(projection=ss_br)

    # Plot the source surface map
    ss_br.plot()
    # Plot the polarity inversion line
    ax.plot_coord(pfss_out.source_surface_pils[0])
    # Plot formatting
    plt.colorbar()
    ax.set_title('Source surface magnetic field')

    plt.show()
    ns = 15
    nphi = 360
    nr = 10
    rss - 2.5

    magnetic_field = pfss_out.bg
    b_r = magnetic_field[:, :, 0, 0]
    b_r = np.copy(b_r)
    # plt.figure()
    # plt.imshow(b_r)

    coors_set = []
    # magnetic_fields = pfss_out.get_bvec(coors_set)
    coor = hmi_map.pixel_to_world(10*u.pix, 20*u.pix)
    plt.figure()
    plt.imshow(hmi_map.data)
