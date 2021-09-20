#!/usr/bin/python

from GT3.Core.Functions.ProfileClasses import OneDProfile, BaseMath
from numpy import array, ndarray



class GTHOneDProfile(OneDProfile):

    def __init__(self, rho, a, val, R=None, Z=None, docs="", units="", plotTitle="", xLabel=r"$\rho$", yLabel="",
                 raw=None, *args, **kwargs):
        """

        :param psi: This is the Psi class that will be used to
        :param val:
        :param R:
        :param Z:
        :param docs:
        :param units:
        :param plotTitle:
        :param xLabel:
        :param yLabel:
        :param wall:
        :param raw:
        """
        super(OneDProfile, self).__init__()
        super(BaseMath, self).__init__()
        __array_priority__ = 1000
        self.__array_priority__ = 1000

        self._data_overwritten = False
        # If there are no values given, we set this to 0s
        # if not np.any(val):
        #     val = np.zeros(rho)
        self.val = val
        self._R = R
        self._Z = Z
        self.a = a

        if type(rho) is list:
            self._rho1D = array(rho)
        elif type(rho) is ndarray:
            self._rho1D = rho

        if type(val) is list:
            self._rho1D = array(val)
        elif type(val) is ndarray:
            self.val = val
        elif type(val) is OneDProfile:
            self.val = val.val

        self._docs = docs
        self._psi = None
        self.units = units
        self.xLabel, self.yLabel, self.plotTitle = xLabel, yLabel, plotTitle
        self.set_plot_rho1d(self._rho1D)
        self._Spline = None
        self._spline_k = 3
        self._spline_s = 2
