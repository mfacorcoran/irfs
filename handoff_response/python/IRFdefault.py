# @file IRFdefault.py
# @brief define default setup options
#
# $Header$
import os
import numpy as np
#p = os.getcwd().split(os.sep)
#print 'loading setup from %s ' % os.getcwd()
# extract class name from file path
#className = p[len(p) - 1]
#print 'eventClass is %s' % className


class Prune(object):
    fileName = 'goodEvent.root'  # file to create
    branchNames = []
    cuts = ''


class Data(object):
    files = []  # use pruned file in event class all by default
    # these correspond to the three runs at SLAC and UW
    generate_area = 6.0
    generated = []
    logemin = []
    logemax = []

# define default binning as attributes of object Bins


class Bins(object):

    @classmethod
    def set_energy_bins(cls, logemin=None, logemax=None, logedelta=None):
        """Convenience method for initializing energy bin edge vector."""
        if logemin is not None:
            cls.logemin = logemin
        if logemax is not None:
            cls.logemax = logemax
        if logedelta is not None:
            cls.logedelta = logedelta

        cls.energy_bins = int((cls.logemax - cls.logemin) / cls.logedelta)
        cls.energy_bin_edges = (np.arange(cls.energy_bins + 1) *
                                cls.logedelta + cls.logemin).tolist()

    @classmethod
    def set_angle_bins(cls, cthmin=None, cthdelta=None):
        """Convenience method for initializing cosTheta bin edge vector."""
        if cthmin is not None:
            cls.cthmin = cthmin
        if cthdelta is not None:
            cls.cthdelta = cthdelta

        cls.angle_bins = int((1.0 - cls.cthmin) / cls.cthdelta)
        cls.angle_bin_edges = np.arange(
            cls.angle_bins + 1) * cls.cthdelta + cls.cthmin

    logemin = 0.75
    logemax = 6.5
    logedelta = 0.25  # 4 per decade

    cthdelta = 0.1
    cthmin = 0.2

    # no overlap with adjacent bins for energy dispersion fits
    edisp_energy_overlap = 0
    edisp_angle_overlap = 0

    # no overlap with adjacent bins for psf fits
    psf_energy_overlap = 0
    psf_angle_overlap = 0


Bins.set_energy_bins()
Bins.set_angle_bins()


class FisheyeBins(Bins):
    """
    subclass of Bins for finer binning of fisheye correction
    """

    logemin = Bins.logemin
    logemax = Bins.logemax
    logedelta = 0.125
    cthmin = 0.2
    cthdelta = 0.05


FisheyeBins.set_energy_bins()
FisheyeBins.set_angle_bins()


class EffectiveAreaBins(Bins):
    """
    subclass of Bins for finer binning of effective area
    """

    @classmethod
    def set_energy_bins(cls, logemin=None, logemax=None):

        if logemin is not None:
            cls.logemin = logemin
        if logemax is not None:
            cls.logemax = logemax

        cls.energy_bin_edges = []
        edges_lo = cls.logemin + \
            np.arange(cls.energy_bins * cls.ebinfactor + 1) * \
            cls.logedelta / cls.ebinfactor
        edges_lo = edges_lo[edges_lo <= cls.ebreak]
        nbin_hi = int(
            (cls.logemax - edges_lo[-1]) / cls.logedelta) * cls.ebinhigh
        edges_hi = edges_lo[-1] + \
            np.arange(nbin_hi + 1) * cls.logedelta / cls.ebinhigh
        cls.energy_bin_edges = list(np.concatenate((edges_lo, edges_hi[1:])))

    logemin = Bins.logemin
    logemax = Bins.logemax
    logedelta = Bins.logedelta
    ebreak = 4.25
    ebinfactor = 4
    ebinhigh = 2


class PSF(object):
    fit_pars = {"ncore": [0.3, 1e-6, 1.],
                "ntail": [0.01, 1e-6, 1.0],
                "score": [0.44, 0.1, 5.0],
                "stail": [1.4, 0.1, 5],
                "gcore": [3.4, 1.001, 20.],
                "gtail": [2.4, 1.001, 20.]}


class Edisp(object):
    Version = 2
    # Fit Parameters key=name, value=[pinit,pmin,max]
    fit_pars = {"f": [0.8, 0.3, 1.0],
                "s1": [1.5, 0.1, 5.0],
                "k1": [1.0, 0.1, 3.0],
                "bias": [0.0, -3, 3],
                "bias2": [0.0, -3, 3],
                "s2": [4.0, 1.2, 10],
                "k2": [1.0, 0.1, 3.0],
                "pindex1": [2.0, 0.1, 5],
                "pindex2": [2.0, 0.1, 5]}


# the log file - to cout if null
logFile = 'log.txt'
