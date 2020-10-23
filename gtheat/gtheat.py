#!/usr/bin/python

import GT3
from scipy.interpolate import UnivariateSpline
import numpy as np
from GT3.utilities.PlotBase import PlotBase
import GT3.constants as constants
from GT3.utilities.PlotBase import MARKERSIZE

e_d = constants.elementary_charge
m_d = constants.deuteron_mass

class chi_i(PlotBase):

    def __init__(self, inputFile):

        self.shot = GT3.gt3(inputFile=inputFile)
        self.shot.run_radial_transport()
        self.rho = self.shot.rtrans.rhor
        self.T = self.shot.rtrans.core.T_fsa
        self.n = self.shot.rtrans.core.n_fsa
        self.Zeff = self.shot.core.z_eff_fsa
        self.nuij = self.shot.rtrans.nu_c_j_k
        self.nuii = self.shot.rtrans.nu_c_j_j
        self.q = self.shot.rtrans.core.q_1D
        self.R = self.shot.core.R[:, 0]
        self.L = self.shot.core.L_fsa
        self.epsilon = self.shot.rtrans.rhor * self.shot.core.a / self.shot.core.R0_a
        self.vth = np.sqrt(self.shot.rtrans.vpol_D**2 + self.shot.rtrans.vtor_D_total**2)
        self.shav_shift = self.shot.core.shaf_shift * (1 - self.rho**2)
        self.dShav_shift = UnivariateSpline(self.rho, self.shav_shift).derivative()(self.rho)
        self.Bp = self.shot.core.B_p_fsa
        self.Bt = self.shot.core.B_t_fsa
        self.set_plot_rho1d(self.rho)
        self.gyro_rad_ion = 1.5E-4 * np.sqrt(self.T.i.ev)/self.Bt
        self.gyro_rad_pol = self.gyro_rad_ion * (self.Bp / self.Bt)
        self.a = self.shot.core.a

        """Do all calculations"""
        self.itg_crit = self._get_itg_crit()
        self.neo_chi = self.calc_neo_ch()
        self.chi_itg_TS_32 = self.calc_itg_TS_32()
        self.chi_itg_TS_12 = self.calc_itg_TS_12()
        self.chi_itg_simp = self.calc_itg_simp()

    def calc_neo_ch(self):

        self.neo_a1, self.neo_a2 = self._neo_ch_get_a(self.n, self.nuij, self.q, self.R, self.epsilon, self.vth)
        self.neo_g1, self.neo_g2 = self._neo_ch_get_g(self.dShav_shift, self.epsilon)
        return np.sqrt(self.epsilon) * self.gyro_rad_pol**2 * self.nuii * (self.neo_a1 * self.neo_g1 + self.neo_a2 * (self.neo_g1 - self.neo_g2))

    def calc_itg_simp(self):

        self.chi_itg_simp = 1.25 * np.sqrt(1 / (self.R * self.L.T.i)) *\
                       (self.gyro_rad_ion * self.T.e.J / (e_d * self.Bt))\
                       * np.heaviside((self.R / self.L.T.i) - self.itg_crit, 0)

        self.chi_itg_simp[self.chi_itg_simp == 0.] = np.nan

    def calc_itg_TS_32(self):
        """
        Calculation of ITG chi using microturbulence and Tore Supra empirical fitting with inverse
        of maximum growth rate. See Eq. 11.62b in Stacey's text
        """

        Ci = 0.014
        q = self.q
        Te = self.T.e.J
        ps = self.gyro_rad_ion
        LTi = self.L.T.i
        R = self.R
        B = self.Bt

        return Ci * q**2 * (Te / (e_d * B)) * (ps / LTi) * (R / LTi)**1.5 *\
               np.heaviside((R / LTi) - self.itg_crit, 0)


    def calc_itg_TS_12(self):
        """
        Calculation of ITG chi using microturbulence and Tore Supra empirical fitting with inverse
        of linar growth rate. See Eq. 11.62b in Stacey's text
        """

        Ci = 0.014
        q = self.q
        Te = self.T.e.J
        ps = self.gyro_rad_ion
        LTi = self.L.T.i
        R = self.R
        B = self.Bt

        return Ci * q ** 2 * (Te / (e_d * B)) * (ps / LTi) * (R / LTi) **.5*\
               np.heaviside((R / LTi) - self.itg_crit, 0)


    def _get_itg_crit(self):

        tau = self.Zeff * self.T.e.ev / self.T.i.ev
        self.dqdr = UnivariateSpline(self.rho, self.q).derivative()(self.rho)
        max_1 = .8 * self.R / self.L.n.e
        max_2 = (1. + 1. / tau) * (1.33 + self.dqdr * 1.91 * self.rho * self.a / (self.q**2)) * (1 - 1.15 * self.epsilon)
        return np.maximum(max_1, max_2)

    def plot_neo_ch(self, edge=True):
        return self._plot_base(self.neo_chi, yLabel="$\chi^{neo}_{r,i}$", edge=edge)

    def plot_itg_chi_simp(self, edge=True):
        return self._plot_base(self.chi_itg_simp, yLabel=r"$\chi^{itg}_{r,i}$", edge=edge)

    def _neo_ch_get_a(self, n, nuij, q, R, eps, vthD):
        """
        Calculate neoclassical chi alpha parameters
        :param n: The densities
        :type n:  GT3.Core.n
        :param nuij: The nu_i_j
        :param q: The safety factor
        :param R: The plasma major radius
        :param eps: The epsilon
        :param vthD: The deuterium poloidal velocity
        """

        alpha = n.C * 36. / n.i
        mustar_i = np.abs(nuij * q * R / (eps**1.5 * vthD))

        a1_numerator = 0.66 * ( 1 + 1.54 * alpha) + (1.88 * np.sqrt(eps)  - 1.54 * eps) * (1 + 3.75 * alpha)
        a1_denominator = 1 + 1.03 * np.sqrt(mustar_i) + 0.31 * mustar_i

        a2 = (0.59 * mustar_i * eps / (1. + 0.74 * mustar_i * eps**1.5)) * (1. + (1.33 * alpha * (1. + 0.6 * alpha) / (1. + 1.79 * alpha)))
        return a1_numerator / a1_denominator, a2

    def _neo_ch_get_g(self, dShift, eps):

        g1 = (1. + 1.5 * (eps**2 + eps * dShift) + (3./8.) * eps**3 * dShift) / (1. + .5 * eps * dShift)
        g2 = (np.sqrt(1. - eps**2) * (1. + (eps * dShift / 2.))) / (1. + (dShift * (np.sqrt(1 - eps**2) - 1.)) / eps)

        return g1, g2

    def plot_T(self, edge=True):
        fig = self._plot_base(self.T.i, yLabel=r"$T_{i,e}[eV]$", edge=edge)
        fig.scatter(self.rho, self.T.e, s=MARKERSIZE)
        fig.legend([r"$T_i$", r"T_e"])

    def plot_n(self, edge=True):
        fig = self._plot_base(self.n.i, yLabel=r"$n_{i,e}[m^{-3}]$", edge=edge)
        fig.scatter(self.rho, self.n.e, s=MARKERSIZE)
        fig.legend([r"$n_i$", r"n_e"])

    def plot_nu(self, edge=True):
        fig = self._plot_base(self.nuij, yLabel=r"$\nu [s^{-1}]$", edge=edge)
        fig.scatter(self.rho, self.nuii,  s=MARKERSIZE)
        fig.legend([r"$\nu_{i,j}$", r"$\nu_{i,i}$"])
        return fig

    def plot_q(self, edge=True):
        return self._plot_base(self.q, yLabel=r"$q$", edge=edge)

    def plot_dqdr(self, edge=True):
        return self._plot_base(self.dqdr, yLabel=r"$\frac{dq}{dr}$", edge=edge)

    def plot_eps(self, edge=True):
        return self._plot_base(self.epsilon, yLabel=r"$\epsilon$", edge=edge)

    def plot_vthermal_D(self, edge=True):
        return self._plot_base(self.vth, yLabel=r"$V_{\theta}[m/s]$", edge=edge)

    def plot_dshift(self, edge=True):
        return self._plot_base(self.dShav_shift, yLabel=r"$\Delta'$", edge=edge)

    def plot_B(self, edge=True):
        fig = self._plot_base(self.Bt, edge=edge)
        fig.scatter(self.rho, self.Bp,  s=MARKERSIZE)
        fig.legend([r"$B_{\phi}$", r"$B_{\theta}$"])
        return fig

    def plot_pol_gyroR(self, edge=True):
        self._plot_base(self.gyro_rad_pol, yLabel=r"$r_{\theta}$", edge=edge)

    def plot_itg_heaviside(self, edge=True):
        itg_crit = self._get_itg_crit()
        fig = self._plot_base((self.R / self.L.T.i), edge=edge)
        fig.scatter(self.rho, itg_crit)
        return fig

    def plot_chis_itg_simp(self, edge=True):
        neoitg = self.neo_chi + self.chi_itg_simp
        fig = self._plot_base(self.shot.rtrans.chi.i.chi1, edge=edge)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi2, color="blue", s=MARKERSIZE)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi4, color="yellow", s=MARKERSIZE)
        fig.scatter(self.rho, self.neo_chi, color="purple", s=MARKERSIZE)
        fig.scatter(self.rho, self.chi_itg_simp, color="green", s=MARKERSIZE)
        fig.scatter(self.rho, neoitg, color="orange", s=MARKERSIZE)
        fig.legend([r"$\chi^{1}_{r,i}$",
                    r"$\chi^{2}_{r,i}$",
                    r"$\chi^{4}_{r,i}$",
                    r"$\chi^{neo}_{r,i}$",
                    r"$\chi^{itg}_{r,i}$",
                    r"$\chi^{neo + itg}_{r,i}"], fontsize=20)
        return fig

    def plot_chis_itg_TS_12(self, edge=True):
        neoitg = self.neo_chi + self.chi_itg_TS_12
        fig = self._plot_base(self.shot.rtrans.chi.i.chi1, edge=edge)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi2, color="blue", s=MARKERSIZE)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi4, color="yellow", s=MARKERSIZE)
        fig.scatter(self.rho, self.neo_chi, color="purple", s=MARKERSIZE)
        fig.scatter(self.rho, self.chi_itg_TS_12, color="green", s=MARKERSIZE)
        fig.scatter(self.rho, neoitg, color="orange", s=MARKERSIZE)
        fig.legend([r"$\chi^{1}_{r,i}$",
                    r"$\chi^{2}_{r,i}$",
                    r"$\chi^{4}_{r,i}$",
                    r"$\chi^{neo}_{r,i}$",
                    r"$\chi^{itg_{TS-1/2}}_{r,i}$",
                    r"$\chi^{neo + itg_{TS-1/2}}_{r,i}$"], fontsize=16)
        return fig

    def plot_chis_itg_TS_32(self, edge=True):
        neoitg = self.neo_chi + self.chi_itg_TS_32
        fig = self._plot_base(self.shot.rtrans.chi.i.chi1, edge=edge)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi2, color="blue", s=MARKERSIZE)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi4, color="yellow", s=MARKERSIZE)
        fig.scatter(self.rho, self.neo_chi, color="purple", s=MARKERSIZE)
        fig.scatter(self.rho, self.chi_itg_TS_32, color="green", s=MARKERSIZE)
        fig.scatter(self.rho, neoitg, color="orange", s=MARKERSIZE)
        fig.legend([r"$\chi^{1}_{r,i}$",
                    r"$\chi^{2}_{r,i}$",
                    r"$\chi^{4}_{r,i}$",
                    r"$\chi^{neo}_{r,i}$",
                    r"$\chi^{itg_{TS-3/2}}_{r,i}$",
                    r"$\chi^{neo + itg_{TS-3/2}}_{r,i}$"], fontsize=16)
        return fig

if __name__ == "__main__":
    import os
    chi_i = chi_i(os.getcwd() + "/inputs/togt3_d3d_123302_2810")
    chi_i.calc_neo_ch()
    chi_i.calc_itg_simp()
    chi_i.calc_itg_TS_12()
    chi_i.calc_itg_TS_32()
    chi_i.plot_chis_itg_TS_32()
    chi_i.plot_chis_itg_TS_12()

