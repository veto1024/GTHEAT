#!/usr/bin/python

import GT3
from scipy.interpolate import UnivariateSpline
import numpy as np
from GT3.utilities.PlotBase import PlotBase
import GT3.constants as constants
from GT3.utilities.PlotBase import MARKERSIZE

e_d = constants.elementary_charge
m_d = constants.deuteron_mass
m_e = constants.electron_mass
mu_0 = constants.mu_0

class chi_i(PlotBase):

    def __init__(self, shot: GT3.gt3):

        self.shot = shot
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
        """Ion thermal velocity """
        self.shav_shift = self.shot.core.shaf_shift * (1 - self.rho**2)
        self.dShav_shift = UnivariateSpline(self.rho, self.shav_shift).derivative()(self.rho)
        self.Bp = self.shot.core.B_p_fsa
        self.Bt = self.shot.core.B_t_fsa
        self.set_plot_rho1d(self.rho)
        self.gyro_rad_ion = self.vth  * m_d / (e_d * self.Bt)
        """Ion gyro radius"""

        self.gyro_rad_pol = self.gyro_rad_ion * (self.Bp / self.Bt)
        self.a = self.shot.core.a
        self.nu_ie = self.shot.rtrans.nu_c_j_e
        self.nu_ei = self.shot.rtrans.nu_c_e_j
        self.nu_ee = self.shot.rtrans.nu_c_e_e
        self.cs = np.sqrt(self.T.e.J/m_d)
        """Ion sound speed"""

        """Do all calculations"""
        self.itg_crit = self._get_itg_crit()
        self.neo_chi = self.calc_neo_ch()
        self.gyro_bohm = self.calc_gyro_bohm()
        self.chi_itg_TS_32 = self.calc_itg_TS_32()
        self.chi_itg_TS_12 = self.calc_itg_TS_12()
        self.chi_itg_simp = self.calc_itg_simp()
        self.chi_DA = self.calc_DA_chi()


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

    def calc_gyro_bohm(self):

        cs = self.cs
        Lpi = self.L.p.i
        ion_r = self.gyro_rad_ion
        return ion_r**2 * cs / Lpi


    def _get_itg_crit(self):

        tau = self.Zeff * self.T.e.ev / self.T.i.ev
        self.dqdr = UnivariateSpline(self.rho, self.q).derivative()(self.rho)
        max_1 = .8 * self.R / self.L.n.e
        max_2 = (1. + 1. / tau) * (1.33 + self.dqdr * 1.91 * self.rho * self.a / (self.q**2)) * (1 - 1.15 * self.epsilon)
        return np.maximum(max_1, max_2)

    def calc_DA_chi(self):
        """
        Calculates Chi from Drift Alfven modes based on gyro-Bohm thermal conductivity
        (see Eq. 11.78 in Stacey).
        :return: Drift Alfven chi
        """
        q = self.q
        R = self.R
        B = self.Bt
        Lpi = self.L.p.i
        ne = self.n.e
        Te = self.T.e.J
        Ti = self.T.i.J
        lambda_e = self.vth / self.nu_ei
        nu_n = (m_d/ m_e)**.25 * ((q * R * Lpi)**.5/lambda_e)
        beta = 2. * mu_0 * ne * Te / B**2
        beta_n = (m_d / m_e)**.5 * (q * R / Lpi) * beta
        k_parallel = 1. / (q * R)
        mu = -1. * k_parallel * Lpi * np.sqrt((m_d * Te)/(m_e * Ti))

        chi_perp_db = (((1 + beta_n**2)**(-3) + nu_n**2) / (1 + beta_n**2 + nu_n**(4./3.)))**.5
        chi_gb = self.calc_gyro_bohm()

        return chi_gb * chi_perp_db / np.sqrt(abs(mu))



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

    def plot_chi_bohm(self, edge=True, show=True):
        return self._plot_base(self.gyro_bohm, yLabel="$\chi^{GB}_{r,i}$", edge=edge, show=show)

    def plot_DA_chi(self, edge=True, show=True):
        return self._plot_base(self.chi_DA, yLabel="$\chi^{DA}_{r,i}$", edge=edge, show=show)

    def plot_neo_ch(self, edge=True, show=True):
        return self._plot_base(self.neo_chi, yLabel="$\chi^{neo}_{r,i}$", edge=edge, show=show)

    def plot_itg_chi_simp(self, edge=True, show=True):
        return self._plot_base(self.chi_itg_simp, yLabel=r"$\chi^{itg}_{r,i}$", edge=edge, show=show)

    def plot_L_T(self, edge=True, show=True):
        fig = self._plot_base(self.L.T.i, yLabel=r"$\L_{T}[m]$", edge=edge, show=show)
        fig.scatter(self.rho, self.L.T.e, s=MARKERSIZE)
        fig.legend([r"$L_{T,i}$", r"$L_{T,e}$"])
        return fig

    def plot_L_n(self, edge=True, show=True):
        fig = self._plot_base(self.L.n.i, yLabel=r"$\L_{n}[m]$", edge=edge, show=show)
        fig.scatter(self.rho, self.L.n.e, s=MARKERSIZE)
        fig.legend([r"$L_{n,i}$", r"$L_{n,e}$"])
        return fig

    def plot_L_P(self, edge=True, show=True):
        fig = self._plot_base(self.L.p.i, yLabel=r"$\L_{P}[m]$", edge=edge, show=show)
        fig.scatter(self.rho, self.L.p.e, s=MARKERSIZE)
        fig.legend([r"$L_{P,i}$", r"$L_{P,e}$"])
        return fig

    def plot_T(self, conv="eV", edge=True, show=True):
        if conv == "eV":
            fig = self._plot_base(self.T.i.ev, yLabel=r"$T_{i,e}[eV]$", edge=edge, show=show)
            fig.scatter(self.rho, self.T.e.ev, s=MARKERSIZE)

        elif conv == "J":
            fig = self._plot_base(self.T.i.J, yLabel=r"$T_{i,e}[J]$", edge=edge, show=show)
            fig.scatter(self.rho, self.T.e.J, s=MARKERSIZE)

        elif conv == "keV":
            fig = self._plot_base(self.T.i.kev, yLabel=r"$T_{i,e}[keV]$", edge=edge, show=show)
            fig.scatter(self.rho, self.T.e.kev, s=MARKERSIZE)

        else:
            print("Conversion key unrecognized. Use 'eV', 'J', or 'keV'.")
            return None

        fig.legend([r"$T_i$", r"T_e"])
        return fig

    def plot_n(self, edge=True, show=True):
        fig = self._plot_base(self.n.i, yLabel=r"$n_{i,e}[m^{-3}]$", edge=edge, show=show)
        fig.scatter(self.rho, self.n.e, s=MARKERSIZE)
        fig.legend([r"$n_i$", r"n_e"])
        return fig

    def plot_nu(self, edge=True, show=True):
        fig = self._plot_base(self.nuij, yLabel=r"$\nu [s^{-1}]$", edge=edge, show=show)
        fig.scatter(self.rho, self.nuii,  s=MARKERSIZE)
        fig.legend([r"$\nu_{i,j}$", r"$\nu_{i,i}$"])
        return fig

    def plot_q(self, edge=True, show=True):
        return self._plot_base(self.q, yLabel=r"$q$", edge=edge, show=show)

    def plot_dqdr(self, edge=True, show=True):
        return self._plot_base(self.dqdr, yLabel=r"$\frac{dq}{dr}$", edge=edge, show=show)

    def plot_eps(self, edge=True, show=True):
        return self._plot_base(self.epsilon, yLabel=r"$\epsilon$", edge=edge, show=show)

    def plot_vthermal_D(self, edge=True, show=True):
        return self._plot_base(self.vth, yLabel=r"$V_{\theta}[m/s]$", edge=edge, show=show)

    def plot_dshift(self, edge=True, show=True):
        return self._plot_base(self.dShav_shift, yLabel=r"$\Delta'$", edge=edge, show=show)

    def plot_B(self, edge=True, show=True):
        fig = self._plot_base(self.Bt, edge=edge, show=show)
        fig.scatter(self.rho, self.Bp,  s=MARKERSIZE)
        fig.legend([r"$B_{\phi}$", r"$B_{\theta}$"])
        return fig

    def plot_pol_ion_gyroR(self, edge=True, show=True):
        return self._plot_base(self.gyro_rad_pol, yLabel=r"$r_{L, \theta}$", edge=edge, show=show)


    def plot_ion_gyroR(self, edge=True, show=True):
        return self._plot_base(self.gyro_rad_ion, yLabel=r"$r_{L}$", edge=edge, show=show)

    def plot_ion_speed(self, edge=True, show=True):
        return self._plot_base(self.cs, yLabel=r"$cs$", edge=edge, show=show)

    def plot_itg_heaviside(self, edge=True, show=True):
        itg_crit = self._get_itg_crit()
        fig = self._plot_base((self.R / self.L.T.i), edge=edge, show=show)
        fig.scatter(self.rho, itg_crit)
        return fig

    def plot_chis_itg_simp(self, edge=True, show=True):
        neoitg = self.neo_chi + self.chi_itg_simp
        fig = self._plot_base(self.shot.rtrans.chi.i.chi1, edge=edge, show=show)
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

    def plot_chis_itg_TS_12(self, edge=True, show=True):
        neoitg = self.neo_chi + self.chi_itg_TS_12
        fig = self._plot_base(self.shot.rtrans.chi.i.chi1, edge=edge, show=show)
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

    def plot_chis_itg_TS_32(self, edge=True, show=True):
        neoitg = self.neo_chi + self.chi_itg_TS_32
        fig = self._plot_base(self.shot.rtrans.chi.i.chi1, edge=edge, show=show)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi2, color="blue", s=MARKERSIZE)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi4, color="yellow", s=MARKERSIZE)
        fig.scatter(self.rho, self.neo_chi, color="purple", s=MARKERSIZE)
        fig.scatter(self.rho, self.chi_itg_TS_32, color="green", s=MARKERSIZE)
        fig.scatter(self.rho, neoitg, color="orange", s=MARKERSIZE)
        fig.legend([r"$\chi^{1}_{r,i}$",
                    r"$\chi^{2}_{r,i}$",
                    r"$\chi^{Stacey}_{r,i}$",
                    r"$\chi^{neo}_{r,i}$",
                    r"$\chi^{itg_{TS-3/2}}_{r,i}$",
                    r"$\chi^{neo + itg_{TS-3/2}}_{r,i}$"], fontsize=16)
        return fig

    def plot_chis_sep(self, edge=True, show=True):
        fig = self._plot_base(self.shot.rtrans.chi.i.chi2, edge=edge, show=show)
        fig.scatter(self.rho, self.shot.rtrans.chi.i.chi4, color="yellow", s=MARKERSIZE)
        fig.scatter(self.rho, self.neo_chi, color="purple", s=MARKERSIZE)
        fig.scatter(self.rho, self.chi_itg_TS_12, color="green", s=MARKERSIZE)
        fig.scatter(self.rho, self.chi_DA, color="orange", s=MARKERSIZE)
        fig.legend([r"$\chi^{2}_{r,i}$",
                    r"$\chi^{Stacey}_{r,i}$",
                    r"$\chi^{neo}_{r,i}$",
                    r"$\chi^{itg_{TS-1/2}}_{r,i}$",
                    r"$\chi^{DA}_{r,i}$"], fontsize=16)
        return fig

    def plot_chis_custom(self, edge=True, neo=True, bohm=False, itg12=False, itg32=False, DA=False, sum=False, show=True):
        """
        Plots customized Chi figure
        :param edge: Where the edge is plotted. Default: True
        :param neo: Whether the neoclassical chi is plotted. Default: False
        :param bohm: Whether the gyro-bohm chi is plotted. Default: False
        :param itg12: Whether the ITG(3/2) chi is plotted. Default: False
        :param itg32: Whether the ITG(3/2) chi is plotted. Default: False
        :param DA: Whether the drift-alfven chi is plotted. Default: False
        :return: The plotted Figure
        """
        if not show:
            yvals = np.array([0.] * len(self.rho))
            if sum:
                if neo:
                    yvals += self.neo_chi
                if bohm:
                    yvals += self.gyro_bohm
                if itg12:
                    yvals += self.chi_itg_TS_12
                if itg32:
                    yvals += self.chi_itg_TS_32
                if DA:
                    yvals += self.chi_DA
                return yvals
            else:
                if neo:
                    return self.neo_chi
                if bohm:
                    return self.gyro_bohm
                if itg12:
                    return self.chi_itg_TS_12
                if itg32:
                    return self.chi_itg_TS_32
                if DA:
                    return self.chi_DA
        else:
            legend = [r"$\chi^{2}_{r,i}$", r"$\chi^{Stacey}_{r,i}$"]
            fig = self._plot_base(self.shot.rtrans.chi.i.chi2, edge=edge, show=show)
            fig.scatter(self.rho, self.shot.rtrans.chi.i.chi4, color="yellow", s=MARKERSIZE)

            if neo:
                fig.scatter(self.rho, self.neo_chi, color="purple", s=MARKERSIZE)
                legend.append(r"$\chi^{neo}_{r,i}$")

            if bohm:
                fig.scatter(self.rho, self.gyro_bohm, color="black", s=MARKERSIZE)
                legend.append(r"$\chi^{GB}_{r,i}$")

            if itg12:
                fig.scatter(self.rho, self.chi_itg_TS_12, color="green", s=MARKERSIZE)
                legend.append(r"$\chi^{itg_{TS-1/2}}_{r,i}$")

            if itg32:
                fig.scatter(self.rho, self.chi_itg_TS_32, color="blue", s=MARKERSIZE)
                legend.append(r"$\chi^{itg_{TS-3/2}}_{r,i}$")

            if DA:
                fig.scatter(self.rho, self.chi_DA, color="orange", s=MARKERSIZE)
                legend.append(r"$\chi^{DA}_{r,i}$")

            fig.legend(legend, fontsize=16)

            return fig