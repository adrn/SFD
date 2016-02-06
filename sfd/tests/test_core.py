""" ...explain... """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os
import sys

# Third-party
from astropy import log as logger
import astropy.coordinates as coord
import astropy.units as u
import numpy as np

from ..core import sfd_ebv, sfd_reddening

def test_sfd_ebv():
    c = coord.SkyCoord(ra=np.random.uniform(0,360,128)*u.degree,
                       dec=np.random.uniform(-10,10,128)*u.degree)
    ebv = sfd_ebv(c)
    assert ebv.shape == (128,)

def test_sfd_reddening_ps1():
    c = coord.SkyCoord(ra=np.random.uniform(0,360,128)*u.degree,
                       dec=np.random.uniform(-10,10,128)*u.degree)
    red_gri = sfd_reddening(c, survey='PS1', filters='gri')
    assert red_gri.shape == (128,3)
