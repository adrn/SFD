""" Note: this was adapted from code by Branimir Sesar """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import json

# Third-party
from astropy import wcs
import astropy.coordinates as coord
from astropy.io import fits
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
import numpy as np
from scipy.ndimage import map_coordinates

__all__ = ['ebv', 'reddening']

def ebv(coordinate, order=1):
    """
    Return SFD E(B-V) at the input coordinate(s).

    Parameters
    ----------
    coordinate : :class:`~astropy.coordinates.SkyCoord`, :class:`~astropy.coordinates.BaseCoordinateFrame`
        The coordinate(s) to compute extinction at.
    order : int (optional)
        Passed to :func:`~scipy.ndimage.map_coordinates`.

    Returns
    -------
    EBV : :class:`~numpy.ndarray`
        Extinction at each input coordinate.

    Example
    -------
    TODO
    h, w = 1000, 4000
    b, l = numpy.mgrid[0:h,0:w]
    l = 180.-(l+0.5) / float(w) * 360.
    b = 90. - (b+0.5) / float(h) * 180.
    ebv = dust.getval(l, b)
    imshow(ebv, aspect='auto', norm=matplotlib.colors.LogNorm())
    """

    ngp_filename = get_pkg_data_filename("data/SFD_dust_4096_ngp.fits")
    sgp_filename = get_pkg_data_filename("data/SFD_dust_4096_sgp.fits")

    # convert input coordinate to Galactic frame
    gal = coordinate.transform_to(coord.Galactic)
    l = gal.l.wrap_at(180*u.degree).degree
    b = gal.b.degree

    # output extinctions
    EBV = np.zeros_like(l, dtype='f4') + np.nan

    b_gtr_zero = b >= 0.
    b_les_zero = np.logical_not(b_gtr_zero)
    if np.any(b_gtr_zero): # north
        hdulist = fits.open(ngp_filename)

        w = wcs.WCS(hdulist[0].header)
        x, y = w.wcs_world2pix(l[b_gtr_zero], b[b_gtr_zero], 0)
        EBV[b_gtr_zero] = map_coordinates(hdulist[0].data, [y, x], order=order, mode='nearest')
        hdulist.close()

    if np.any(b_les_zero): # south
        hdulist = fits.open(sgp_filename)

        w = wcs.WCS(hdulist[0].header)
        x, y = w.wcs_world2pix(l[b_les_zero], b[b_les_zero], 0)
        EBV[b_les_zero] = map_coordinates(hdulist[0].data, [y, x], order=order, mode='nearest')
        hdulist.close()

    return EBV

def reddening(coordinate, survey, filters, order=1):
    """
    Return SFD reddening in the given filters

    Parameters
    ----------
    coordinate : :class:`~astropy.coordinates.SkyCoord`, :class:`~astropy.coordinates.BaseCoordinateFrame`
        The coordinate(s) to compute extinction at.
    survey : str
        Survey name in the E(B-V) filter conversion file.
    filters : iterable, str
    order : int (optional)
        Passed to :func:`~scipy.ndimage.map_coordinates`.

    Returns
    -------
    red : :class:`~numpy.ndarray`
        Reddening in each filter for each input coordinate. Output shape will
        be `(len(coordinate), len(filters))`.

    Example
    -------
    TODO:
    """

    filters = list(filters)

    ebv_conv_file = get_pkg_data_filename("data/EBV_to_filter.json")
    with open(ebv_conv_file,'r') as f:
        ebv_to_red = json.loads(f.read())[survey]

    EBV = sfd_ebv(coordinate, order=order)

    out = np.zeros((len(coordinate), len(filters)))
    for i,filter_name in enumerate(filters):
        out[:,i] = ebv_to_red[filter_name]*EBV

    return out
