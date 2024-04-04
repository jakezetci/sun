import pfsspy
import pfsspy.utils

import os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import frame_transform_graph
import sunpy.coordinates as spc
import astropy.coordinates as coord


class ASineHGC(spc.HeliographicCarrington):
    """
    This frame is a reverse transformation from the HGC frame with CUNIT2 sin(deg) to the real HGC frame.
    It is defined backwards because the WCS object for the HGC sin(deg) frame reads in as a proper HGC frame
    due to the values of CTYPE. This frame therefore converts from what WCS thinks is HGC to ASineHGC but is
    in reality converting from ASineHGC to HGC.
    """


@frame_transform_graph.transform(coord.FunctionTransform, spc.HeliographicCarrington, ASineHGC)
def hgc_to_sinehgc(hgc_coord, sinehgc_frame):
    lat = hgc_coord.lat
    lon = hgc_coord.lon

    lat_out = u.Quantity(np.arcsin(lat.value), u.rad)

    return sinehgc_frame.realize_frame(coord.UnitSphericalRepresentation(lat=lat_out, lon=lon))


@frame_transform_graph.transform(coord.FunctionTransform, ASineHGC, spc.HeliographicCarrington)
def sinehgc_to_hgc(sinehgc_coord, hgc_frame):
    lat = sinehgc_coord.lat
    lon = sinehgc_coord.lon

    lat_out = u.Quantity(np.sin(lat.value), u.deg)

    return hgc_frame.realize_frame(coord.UnitSphericalRepresentation(lat=lat_out, lon=lon))


if __name__ == "__main__":
    date = '2020-02-13'
    date = np.datetime64(date)
    date2 = date + np.timedelta64(1, 'D')
    time = a.Time(date, date2)
    # time = a.Time('2010/01/01', '2010/01/02')
    series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')
    series = a.jsoc.Series('hmi.mldailysynframe_720s')

    # crot = a.jsoc.PrimeKey('CAR_ROT',
    #                       int(sunpy.coordinates.sun.carrington_rotation_number(time)))
    # print(sunpy.coordinates.sun.carrington_rotation_number(time))
    # result = Fido.search(series, crot,
    #                     a.jsoc.Notify('rrzhdanov@edu.hse.ru'))
    result = Fido.search(time, series,
                         a.jsoc.Notify('rrzhdanov@edu.hse.ru'))
    files = Fido.fetch(result)
    hdus = fits.open(files[0])

    data = hdus[0].data
    header = dict(hdus[0].header)
    header['CUNIT2'] = 'deg'
    # header['CDELT2'] = 180 / np.pi * header['CDELT2']
    del header['HGLN_OBS']

    print(files)

    hmi_map = sunpy.map.Map((data, header))
    hmi_map.plot_settings['cmap'] = 'hmimag'
    hmi_map.plot_settings['norm'] = plt.Normalize(-1500, 1500)
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

    print(hmi_map.coordinate_frame)
    # Create a figure with the Map's projection:
    fig = plt.figure(figsize=(12, 5))
    axes = plt.subplot(projection=hmi_map)

    # Plot the image
    im = hmi_map.plot()

    # Set up the Sine Latitude Grid
    x = axes.coords[0]
    y = axes.coords[1]

    x.set_coord_type('longitude', coord_wrap=360.)

    x.set_major_formatter('dd')
    y.set_major_formatter('d.d')

    x.set_axislabel("Carrington Longitude [deg]")
    y.set_axislabel("Sine Latitude")

    x.set_ticks(color='black', exclude_overlapping=True)
    y.set_ticks(color='black', exclude_overlapping=True)

    # Hide the grid
    axes.coords.grid(alpha=0)

    # Create a colorbar
    cb = plt.colorbar(im, fraction=0.019, pad=0.1)
    cb.set_label("LOS Magnetic Field [Gauss]")

    # Now create the overlay in the actual HGC coordinate frame

    overlay = axes.get_coords_overlay('asinehgc')

    lon = overlay[0]
    lat = overlay[1]

    lon.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)

    lat.set_major_formatter('dd')

    lat.set_axislabel('Solar Latitude [deg]')

    lat.set_ticks_position('tr')
    lat.set_ticks(spacing=10*u.deg, exclude_overlapping=True)

    # Another horrible hack to make the ticks draw on the RHS
    axes.set_xlim((0, 3585))

    plt.title("HMI Daily Synoptic Frame for Carrington Rotation {}-{}".format(
        header['CAR_ROT'], header['CAR_ROT']+1))
    """
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
    """
