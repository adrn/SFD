# SFD

Schlegel, Finkbeiner, David dust maps.

## Installation

The dust map data is stored in this repository using [git-lfs](https://git-lfs.github.com). Clone the repository and install using

    python setup.py install

## Dependencies

This requires scipy, astropy, and numpy. Install with

    conda install numpy, scipy, astropy

or

    pip install -r pip-requirements.txt

from the top-level of the cloned repository

## Usage

    import astropy.units as u
    import astropy.coordinates as coord
    from sfd import ebv, reddening

    c = coord.SkyCoord(ra=[154.12, 11.1]*u.degree,
                       dec=[-21.63,31.65]*u.degree)

    EBV = ebv(c)

    # or

    reddening = reddening(c, survey='PS1', filters='gri')

