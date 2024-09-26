from functools import lru_cache
from typing import Tuple
from scipy.optimize import brentq, differential_evolution
import numpy as np
import pandas as pd
import s3fs


class PCSAFT:

    def __init__(self):
        self.param = None
        self.param_init()

    def param_init(self):
        self.param = {
            'N_av': 6.022140857e23,
            'kb': 1.380648465952442093e-23,
            'pi': 3.141592653589793,
            'a0': np.array(
                [0.910563145, 0.636128145, 2.686134789, -26.54736249, 97.75920878, -159.5915409, 91.29777408]),
            'a1': np.array(
                [-0.308401692, 0.186053116, -2.503004726, 21.41979363, -65.25588533, 83.31868048, -33.74692293]),
            'a2': np.array(
                [-0.090614835, 0.452784281, 0.596270073, -1.724182913, -4.130211253, 13.77663187, -8.672847037]),
            'b0': np.array(
                [0.724094694, 2.238279186, -4.002584949, -21.00357682, 26.85564136, 206.5513384, -355.6023561]),
            'b1': np.array(
                [-0.575549808, 0.699509552, 3.892567339, -17.21547165, 192.6722645, -161.8264617, -165.2076935]),
            'b2': np.array(
                [0.097688312, -0.255757498, -9.155856153, 20.64207597, -38.80443005, 93.62677408, -29.66690559]),
            'a0dip': [0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308],
            'a1dip': [0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135],
            'a2dip': [-1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575],
            'b0dip': [0.2187939, -1.1896431, 1.1626889, 0, 0],
            'b1dip': [-0.5873164, 1.2489132, -0.5085280, 0, 0],
            'b2dip': [3.4869576, -14.915974, 15.372022, 0, 0],
            'c0dip': [-0.0646774, 0.1975882, -0.8087562, 0.6902849, 0],
            'c1dip': [-0.9520876, 2.9924258, -2.3802636, -0.2701261, 0],
            'c2dip': [-0.6260979, 1.2924686, 1.6542783, -3.4396744, 0],
            'conv': 7242.702976750923,
            'atomic constants': np.array(
                [1.990, 1.699, -0.205, -0.017, 3.706, -1.693, 8.88, 1.108, 1.548, 0.114, 2.016, 15.226,
                 3.504, 4.83, -2.755, 0.954])
        }

    def compute_z(self, rho: float, T: float, comp: np.ndarray, param):
        d = param.s * (1 - 0.12 * np.exp(-3 * param.e / T))

        Den = rho * self.param['N_av'] / 1e30
        zeta = self.param['pi'] / 6 * Den * np.sum((comp * param.m)[:, np.newaxis] * d[:, np.newaxis] ** np.arange(4),
                                                   axis=0)

        m_avg = np.sum(comp * param.m)

        g_hs, delta_ghs_rho = self._compute_hs_terms(d, zeta)

        s_ij = 0.5 * (param.s[:, np.newaxis] + param.s)
        e_ij = np.sqrt(np.outer(param.e, param.e)) * (1 - param.k_ij)

        m2es3 = np.sum(np.outer(comp * param.m, comp * param.m) * (e_ij / T) * s_ij ** 3)
        m2e2s3 = np.sum(np.outer(comp * param.m, comp * param.m) * (e_ij / T) ** 2 * s_ij ** 3)

        Zhs = self._compute_zhs(zeta)
        sum1 = np.sum(comp * (param.m - 1) / np.diag(g_hs) * np.diag(delta_ghs_rho))

        a, b = self._compute_ab(m_avg)
        I2, deltaI1_det, deltaI2_det = self._compute_I_terms(zeta[3], a, b)

        C1, C2 = self._compute_C_terms(m_avg, zeta[3])

        Zhc = m_avg * Zhs - sum1
        Zdisp = -2 * self.param['pi'] * Den * deltaI1_det * m2es3 - self.param['pi'] * Den * m_avg * (
                C1 * deltaI2_det + C2 * zeta[3] * I2) * m2e2s3

        return 1 + Zhc + Zdisp

    def _compute_hs_terms(self, d: np.ndarray, zeta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        d_ij = d[:, np.newaxis] + d
        d_ij_term = d[:, np.newaxis] * d / d_ij
        g_hs = 1 / (1 - zeta[3]) + d_ij_term * 3 * zeta[2] / (1 - zeta[3]) ** 2 + d_ij_term ** 2 * 2 * zeta[2] ** 2 / (
                1 - zeta[3]) ** 3
        delta_ghs_rho = zeta[3] / (1 - zeta[3]) ** 2 + d_ij_term * (
                3 * zeta[2] / (1 - zeta[3]) ** 2 + 6 * zeta[2] * zeta[3] / (1 - zeta[3]) ** 3) + d_ij_term ** 2 * (
                                4 * zeta[2] ** 2 / (1 - zeta[3]) ** 3 + 6 * zeta[2] ** 2 * zeta[3] / (
                                1 - zeta[3]) ** 4)

        return g_hs, delta_ghs_rho

    def _compute_zhs(self, zeta: np.ndarray):
        return zeta[3] / (1 - zeta[3]) + 3 * zeta[1] * zeta[2] / zeta[0] / (1 - zeta[3]) ** 2 + \
            (3 * (zeta[2] ** 3) - zeta[3] * (zeta[2] ** 3)) / zeta[0] / (1 - zeta[3]) ** 3

    def _compute_ab(self, m_avg: float) -> tuple[np.ndarray, np.ndarray]:
        a = self.param['a0'] + (m_avg - 1) / m_avg * self.param['a1'] + (m_avg - 1) / m_avg * (m_avg - 2) / m_avg * \
            self.param['a2']
        b = self.param['b0'] + (m_avg - 1) / m_avg * self.param['b1'] + (m_avg - 1) / m_avg * (m_avg - 2) / m_avg * \
            self.param['b2']

        return a, b

    def _compute_I_terms(self, eta3: float, a: np.ndarray, b: np.ndarray):
        eta3_powers = eta3 ** np.arange(7)
        I2 = np.sum(b * eta3_powers)
        deltaI1_det = np.sum(a * np.arange(1, 8) * eta3_powers)
        deltaI2_det = np.sum(b * np.arange(1, 8) * eta3_powers)
        return I2, deltaI1_det, deltaI2_det

    def _compute_C_terms(self, m_avg: float, eta3: float) -> tuple[float, float]:
        C1 = 1 / (1 + m_avg * (8 * eta3 - 2 * eta3 ** 2) / ((1 - eta3) ** 4) + (1 - m_avg) * (
                20 * eta3 - 27 * eta3 ** 2 + 12 * eta3 ** 3 - 2 * eta3 ** 4) / (((1 - eta3) * (2 - eta3)) ** 2))
        C2 = -C1 ** 2 * (m_avg * (-4 * eta3 ** 2 + 20 * eta3 + 8) / ((1 - eta3) ** 5) + (1 - m_avg) * (
                2 * eta3 ** 3 + 12 * eta3 ** 2 - 48 * eta3 + 40) / (((1 - eta3) * (2 - eta3)) ** 3))
        return C1, C2

    def compute_density(self, T, p, comp, phase, param):

        try:
            self._validate_inputs(T, p, comp, phase, param)

            num_pts = 50
            rho_guesses = np.linspace(1e-13, 0.7405, num_pts)
            x_low, x_up = self._find_root_intervals(rho_guesses, T, p, comp, param)

            if len(x_low) == 1:
                return self._solve_density(x_low[0], x_up[0], T, p, comp, param)
            elif 1 < len(x_low) <= 3:
                return self._solve_density(x_low[-1] if phase == 0 else x_low[0],
                                           x_up[-1] if phase == 0 else x_up[0],
                                           T, p, comp, param)
            elif len(x_low) > 3:
                return self._solve_multiple_roots(x_low, x_up, T, p, comp, param)
            else:
                return self._global_optimization(T, p, comp, param)
        except Exception as e:
            raise ValueError(f"无法计算密度: {str(e)}")

    def _solve_density(self, rho_low, rho_up, T, p, comp: np.ndarray,
                       param):

        rho_low = self.compute_molar_density(rho_low, T, len(comp), comp, param)
        rho_up = self.compute_molar_density(rho_up, T, len(comp), comp, param)

        def f(rho):
            return self.solve_rho(rho, T, p, comp, param)

        return brentq(f, rho_low, rho_up, xtol=1e-8, maxiter=200)

    def _solve_multiple_roots(self, x_low, x_up, T: float, p: float, comp: np.ndarray,
                              param) -> float:

        bounds = list(zip(x_low, x_up))
        result = differential_evolution(lambda x: self.compute_Gres(T, x[0], comp, param),
                                        bounds, popsize=20, tol=1e-8)
        return result.x[0]

    def _global_optimization(self, T: float, p: float, comp, param) -> float:
        bounds = [(1e-8, 7.4)]
        result = differential_evolution(
            lambda x: abs(self.compute_pressure_residual(x[0], T, p, comp, param)),
            bounds, popsize=20, tol=1e-8)
        return self.compute_molar_density(result.x[0], T, len(comp), comp, param)

    def _find_root_intervals(self, rho_guesses, T: float, p: float, comp: np.ndarray, param):

        x_low, x_up = [], []
        P_err_prev = self.compute_pressure_residual(
            self.compute_molar_density(rho_guesses[0], T, len(comp), comp, param), T, p, comp, param)

        for rho_guess in rho_guesses[1:]:
            rho_molar = self.compute_molar_density(rho_guess, T, len(comp), comp, param)
            P_err = self.compute_pressure_residual(rho_molar, T, p, comp, param)

            if P_err_prev * P_err < 0:
                x_low.append(rho_guesses[rho_guesses < rho_guess][-1])
                x_up.append(rho_guess)

            P_err_prev = P_err

        return x_low, x_up

    def _validate_inputs(self, T: float, p: float, comp: np.ndarray, phase: int, param):
        """验证输入参数"""
        if T <= 0 or p <= 0:
            raise ValueError("温度和压力必须为正值")
        if not (0 <= phase <= 1):
            raise ValueError("相态必须为0（液相）或1（气相）")
        if not (0.99 <= np.sum(comp) <= 1.01):
            raise ValueError("组分之和必须接近1")

    def compute_molar_density(self, rho_low, T, n, comp, param):

        d = param.s * (1 - 0.12 * np.exp(-3 * param.e / T))
        sum_c = np.sum(comp * param.m * d ** 3)
        Molar_Density = 6 / self.param['pi'] * rho_low / sum_c * 1.0e30 / self.param['N_av']
        return Molar_Density

    def compute_pressure_residual(self, rho_guess, T, p, comp, param):

        Pcal = self._compute_P(rho_guess, T, comp, param)
        resid = (Pcal - p) / p
        if np.isinf(resid):
            P_Residual = 10e300
        else:
            P_Residual = resid
        return P_Residual

    def _compute_P(self, rho_guess, T, comp, param):

        Den = rho_guess * self.param['N_av'] / 1e30
        Z = self.compute_z(rho_guess, T, comp, param)
        # pa
        P = Z * self.param['kb'] * T * Den * 1e30

        return P

    def solve_rho(self, rho, T, p, comp, param):

        pcal = self._compute_P(rho, T, comp, param)
        return pcal - p

    def compute_Gres(self, T, rho, comp, param):

        ares = self.compute_Ares(T, rho, comp, param)

        Z = self.compute_z(rho, T, comp, param)

        P = self._compute_P(rho, T, comp, param)

        gres = (ares + (Z - 1) - np.log(Z)) * self.param['kb'] * self.param['N_av'] * T + self.param['kb'] * self.param[
            'N_AV'] * T * np.log(P / 101325)

        Gres = gres

        return Gres

    def compute_Ares(self, T, rho, comp, param):

        d = param.s * (1 - 0.12 * np.exp(-3 * param.e / T))
        Den = rho * self.param['N_av'] / 1e30
        zeta = self.param['pi'] / 6 * Den * np.sum((comp * param.m)[:, np.newaxis] * d[:, np.newaxis] ** np.arange(4),
                                                   axis=0)

        m_avg = np.sum(comp * param.m)

        g_hs, delta_ghs_rho = self._compute_hs_terms(d, zeta)
        s_ij = 0.5 * (param.s[:, np.newaxis] + param.s)
        e_ij = np.sqrt(np.outer(param.e, param.e)) * (1 - param.k_ij)

        m2es3 = np.sum(np.outer(comp * param.m, comp * param.m) * (e_ij / T) * s_ij ** 3)
        m2e2s3 = np.sum(np.outer(comp * param.m, comp * param.m) * (e_ij / T) ** 2 * s_ij ** 3)

        ares_hs = self._compute_Ares_hs(zeta)

        sum3 = np.sum(comp * (param.m - 1) * np.log(np.diag(g_hs)))

        a, b = self._compute_ab(m_avg)

        eta3_powers = zeta[3] ** np.arange(7)
        I1 = np.sum(a * eta3_powers)
        I2 = np.sum(b * eta3_powers)
        C1 = 1 / (1 + m_avg * (8 * zeta[3] - 2 * zeta[3] ** 2) / ((1 - zeta[3]) ** 4) + (1 - m_avg) * (
                20 * zeta[3] - 27 * zeta[3] ** 2 + 12 * (zeta[3] ** 3) - 2 * (zeta[3] ** 4)) / (
                      (((1 - zeta[3]) * (2 - zeta[3])) ** 2)))

        ares_hc = m_avg * ares_hs - sum3
        ares_disp = -2 * self.param['pi'] * Den * I1 * m2es3 - self.param['pi'] * Den * m_avg * C1 * I2 * m2e2s3;
        ares = ares_hc + ares_disp
        return ares

    def _compute_Ares_hs(self, zeta):

        return 1 / zeta[0] * (
                3 * zeta[1] * zeta[2] / (1 - zeta[3]) + (zeta[2] ** 3) / (zeta[3] * (1 - zeta[3]) ** 2) + (
                (zeta[2] ** 3) / (zeta[3] ** 2) - zeta[0]) * np.log(1 - zeta[3]))

    def compute_hres(self, T, rho, comp, param):

        Z = self.compute_z(rho, T, comp, param)
        dares_dt = self.compute_dadT(T, rho, comp, param)
        hres = (-T * dares_dt + (Z - 1)) * self.param['kb'] * self.param['N_av'] * T

        return hres

    def compute_dadT(self, T, rho, comp, param):

        d = param.s * (1 - 0.12 * np.exp(-3 * param.e / T))
        dd_dT = param.s * -3 * param.e / T ** 2 * 0.12 * np.exp(-3 * param.e / T)
        den = rho * self.param['N_av'] / 1e30
        zeta = self.param['pi'] / 6 * den * np.sum((comp * param.m)[:, np.newaxis] * d[:, np.newaxis] ** np.arange(4),
                                                   axis=0)
        eta3 = zeta[3]
        m_avg = np.sum(comp * param.m)
        dzeta_dT = self.param['pi'] / 6 * den * np.sum(
            (comp * param.m)[:, np.newaxis] * np.arange(4)[np.newaxis, :] * dd_dT[:, np.newaxis] * d[:, np.newaxis] ** (
                    np.arange(4) - 1), axis=0)

        ghs, delta_ghs_rho = self._compute_hs_terms(d, zeta)
        ddij_dt = (d[:, np.newaxis] * d[np.newaxis, :] / (d[:, np.newaxis] + d[np.newaxis, :])) * (
                dd_dT[:, np.newaxis] / d[:, np.newaxis] + dd_dT[np.newaxis, :] / d[np.newaxis, :] -
                (dd_dT[:, np.newaxis] + dd_dT[np.newaxis, :]) / (d[:, np.newaxis] + d[np.newaxis, :]))
        s_ij = 0.5 * (param.s[:, np.newaxis] + param.s[np.newaxis, :])
        e_ij = np.sqrt(param.e[:, np.newaxis] * param.e[np.newaxis, :]) * (1 - param.k_ij)

        m2es3 = np.sum(comp[:, np.newaxis] * comp[np.newaxis, :] * param.m[:, np.newaxis] * param.m[np.newaxis, :] * (
                e_ij / T) * s_ij ** 3)
        m2e2s3 = np.sum(comp[:, np.newaxis] * comp[np.newaxis, :] * param.m[:, np.newaxis] * param.m[np.newaxis, :] * (
                (e_ij / T) ** 2) * s_ij ** 3)
        dghs_dt = (dzeta_dT[3] / (1 - zeta[3]) ** 2 +
                   3 * (ddij_dt * zeta[2] + (
                        d[:, np.newaxis] * d[np.newaxis, :] / (d[:, np.newaxis] + d[np.newaxis, :])) * dzeta_dT[
                            2]) / (1 - zeta[3]) ** 2 +
                   4 * (d[:, np.newaxis] * d[np.newaxis, :] / (d[:, np.newaxis] + d[np.newaxis, :])) * zeta[2] * (
                           1.5 * dzeta_dT[3] + ddij_dt * zeta[2] +
                           (d[:, np.newaxis] * d[np.newaxis, :] / (d[:, np.newaxis] + d[np.newaxis, :])) * dzeta_dT[
                               2]) / (1 - zeta[3]) ** 3 +
                   6 * ((d[:, np.newaxis] * d[np.newaxis, :] / (d[:, np.newaxis] + d[np.newaxis, :])) * zeta[2]) ** 2 *
                   dzeta_dT[3] / (1 - zeta[3]) ** 4)

        dadt_hs = (1 / zeta[0] * (3 * (dzeta_dT[1] * zeta[2] + zeta[1] * dzeta_dT[2]) / (1 - zeta[3]) +
                                  3 * zeta[1] * zeta[2] * dzeta_dT[3] / (1 - zeta[3]) ** 2 +
                                  3 * zeta[2] ** 2 * dzeta_dT[2] / zeta[3] / (1 - zeta[3]) ** 2 +
                                  zeta[2] ** 3 * dzeta_dT[3] * (3 * zeta[3] - 1) / zeta[3] ** 2 / (1 - zeta[3]) ** 3 +
                                  (3 * zeta[2] ** 2 * dzeta_dT[2] * zeta[3] - 2 * zeta[2] ** 3 * dzeta_dT[3]) / zeta[
                                      3] ** 3 * np.log(1 - zeta[3]) +
                                  (zeta[0] - zeta[2] ** 3 / zeta[3] ** 2) * dzeta_dT[3] / (1 - zeta[3])))

        a, b = self._compute_ab(m_avg)

        I1 = np.sum(a * eta3 ** np.arange(7))
        I2 = np.sum(b * eta3 ** np.arange(7))
        detI1_det = np.sum(a * np.arange(7) * eta3 ** (np.arange(7) - 1) * dzeta_dT[3])
        detI2_det = np.sum(b * np.arange(7) * eta3 ** (np.arange(7) - 1) * dzeta_dT[3])

        C1 = 1 / (1 + m_avg * (8 * eta3 - 2 * eta3 ** 2) / (1 - eta3) ** 4 +
                  (1 - m_avg) * (20 * eta3 - 27 * eta3 ** 2 + 12 * eta3 ** 3 - 2 * eta3 ** 4) / (
                          (1 - eta3) * (2 - eta3)) ** 2)
        C2 = -C1 ** 2 * (m_avg * (-4 * eta3 ** 2 + 20 * eta3 + 8) / (1 - eta3) ** 5 +
                         (1 - m_avg) * (2 * eta3 ** 3 + 12 * eta3 ** 2 - 48 * eta3 + 40) / (
                                 (1 - eta3) * (2 - eta3)) ** 3)
        dC1_dt = C2 * dzeta_dT[3]

        summ = np.sum(comp * (param.m - 1) * np.diag(dghs_dt) / np.diag(ghs))
        dadt_hc = m_avg * dadt_hs - summ

        dadt_disp = -2 * self.param['pi'] * den * (detI1_det - I1 / T) * m2es3 - self.param['pi'] * den * m_avg * (
                dC1_dt * I2 + C1 * detI2_det - 2 * C1 * I2 / T) * m2e2s3

        return dadt_hs + dadt_hc + dadt_disp

    def compute_sres(self, T, rho, comp, param):

        ares = self.compute_Ares(T, rho, comp, param)
        Z = self.compute_z(rho, T, comp, param)
        P = self._compute_P(rho, T, comp, param)
        Sres = (-T * (self.compute_dadT(T, rho, comp, param) + ares / T) + np.log(
            Z)) * self.param['kb'] * self.param['N_av'] - self.param['kb'] * self.param['N_av'] * np.log(P / 101325)
        return Sres

    def compute_den_correction(self, rho, polymer_zeroth_mole_frac, dpn):

        return rho * (1 - polymer_zeroth_mole_frac) + rho * polymer_zeroth_mole_frac * dpn

    def compute_z_correction(self, z, polymer_first_mole_frac, polymer_zeroth_mole_frac):

        return z * (1 - polymer_first_mole_frac) + z * polymer_zeroth_mole_frac

    def compute_hres_correction(self, hres, dpn):

        return hres / dpn

    def compute_hig(self, param, T):

        flag = param.type

        Tsp = np.linspace(start=298.15, stop=T, num=10001)

        cpig_arr = self.compute_cpig(Tsp[0], param, flag)

        for idx, t in enumerate(Tsp):

            if idx == 0:
                continue
            cpig_arr = np.vstack((cpig_arr, self.compute_cpig(t, param, flag)))

        hig = param.HIGTREF * 1e6 + (T - 298.15) / (2 * 10001) * (
                cpig_arr[0, :] + np.sum(cpig_arr[1:-1, :], axis=0) + cpig_arr[-1, :])

        return hig

    def compute_hig_mix(self, hig_arr, comp):

        return np.sum(hig_arr * comp)

    def compute_cpig(self, T, param, flag):

        cpig = np.zeros(len(flag))
        for i, f in enumerate(flag):

            # 聚合物物质
            if f == 1:

                c7 = param.CPIG[i][6]
                c8 = param.CPIG[i][7]

                if c7 <= T <= c8:
                    cpig[i] = np.sum(T ** np.arange(6) * param.CPIG[i][0:6])
                elif T < c7:
                    cpig[i] = param.CPIG[i][8] + param.CPIG[i][9] * T ** param.CPIG[i][10]
                else:
                    cpig[i] = (np.sum((param.CPIG[i][7] + 1e-6) ** np.arange(6) * param.CPIG[i][0:6]) - np.sum(
                        param.CPIG[i][7] ** np.arange(6) * param.CPIG[i][0:6])) / 1e-6 * (T - c8) + np.sum(
                        param.CPIG[i][7] ** np.arange(6) * param.CPIG[i][0:6])
            elif f == 0:

                cpig[i] = param.CPIG[i][0] + param.CPIG[i][1] * (
                        param.CPIG[i][2] / T / np.sinh(param.CPIG[i][2] / T)) ** 2 + param.CPIG[i][3] * (
                                  param.CPIG[i][4] / T / np.cosh(param.CPIG[i][4] / T)) ** 2

        return cpig

    def compute_Enthalpy(self, hig, hres):

        return hig / 1000 + hres

    def retrive_param_from_DB(self, CAS, param_name):
        fs = s3fs.S3FileSystem()
        db_path = 'db.csv'
        # 需要在IBD新建一个存放文件的文件夹，并将路径进行修改
        db_df = pd.read_csv(fs.open('/data/db.csv'))
        db_df.set_index(db_df.iloc[:, 0], inplace=True)
        return db_df[CAS][param_name]

    def compute_glassify_temperature_Askadskii_Matveev(self, param, co_molefrac):

        if len(param.type) > 2:
            print('暂不支持多元（>2)的共聚体系')
            return 0

        if 1 not in param.type:

            return 0

        else:

            Tg_exp_arr = np.array([], dtype=np.float32)
            van_der_waar_arr = np.array([], dtype=np.float32)
            for idx, name in enumerate(param.CAS):
                Tg_exp_arr = np.append(Tg_exp_arr, np.array(self.retrive_param_from_DB(name, 'Tg'), dtype=np.float32))
                van_der_waar_arr = np.append(van_der_waar_arr,
                                             np.array(self.retrive_param_from_DB(name, 'DVAMVDW'), dtype=np.float32))

            Tg = np.sum(van_der_waar_arr * co_molefrac) / (
                    np.sum(co_molefrac * van_der_waar_arr / Tg_exp_arr) + 0.03 * np.sum(
                co_molefrac * (1 - co_molefrac)))
            return Tg

    def compute_solid_density_Askadskii_Matveev(self, param, T, co_molefrac):

        kg = 0.737
        Tg = self.compute_glassify_temperature_Askadskii_Matveev(param, co_molefrac)
        alpha_L = np.array([], dtype=np.float32)
        alpha_G = np.array([], dtype=np.float32)
        Mw_arr = np.array([], dtype=np.float32)
        van_der_waar_arr = np.array([], dtype=np.float32)
        for idx, i in enumerate(param.CAS):
            alpha_L = np.append(alpha_L, np.array(self.retrive_param_from_DB(i, 'CTE'), dtype=np.float32))
            van_der_waar_arr = np.append(van_der_waar_arr,
                                         np.array(self.retrive_param_from_DB(i, 'DVAMVDW'), dtype=np.float32))
            Mw_arr = np.append(Mw_arr, np.array(self.retrive_param_from_DB(i, 'MW'), dtype=np.float32))
        alpha_G = alpha_L - 0.113 / Tg

        if T < Tg:

            rho = kg * np.sum(Mw_arr * co_molefrac) / (1 + np.sum(alpha_G * co_molefrac) * (T - Tg)) / np.sum(
                van_der_waar_arr * 0.6022 * co_molefrac)
        elif T > Tg:

            rho = kg * np.sum(Mw_arr * co_molefrac) / (1 + np.sum(alpha_L * co_molefrac) * (T - Tg)) / np.sum(
                van_der_waar_arr * 0.6022 * co_molefrac)
        else:
            rho = kg * np.sum(Mw_arr * co_molefrac) / np.sum(
                van_der_waar_arr * 0.6022 * co_molefrac)

        return rho

    def compute_solid_heat_capacity_Bicerano(self, param, T, co_molefrac):

        Tg = self.compute_glassify_temperature_Askadskii_Matveev(param, co_molefrac)

        CPsTREF_exp = np.array([], dtype=np.float32)
        CPlTREF_exp = np.array([], dtype=np.float32)

        for idx, cas in enumerate(param.CAS):
            CPsTREF_exp = np.append(CPsTREF_exp, np.array(self.retrive_param_from_DB(cas, 'CpsTref'), dtype=np.float32))
            CPlTREF_exp = np.append(CPlTREF_exp, np.array(self.retrive_param_from_DB(cas, 'CplTref'), dtype=np.float32))

        if 100 < T <= Tg:
            CP = np.sum(CPsTREF_exp * (1 + 3 * 1e-3 * (T - 298)) * co_molefrac)
        elif T > Tg:
            CP = np.sum(CPlTREF_exp * (1 + 1.3 * 1e-3 * (T - 298)) * co_molefrac)

        return CP

    def compute_solid_thermal_conductivity_Askadskii_Matveev(self, param, co_molefrac, T):

        Tg = self.compute_glassify_temperature_Askadskii_Matveev(param, co_molefrac)
        ac_sum_arr = np.array([], dtype=np.float32)
        vdw_arr = np.array([], dtype=np.float32)
        mw_arr = np.array([], np.float32)
        a_num = np.array([])
        for i, cas in enumerate(param.CAS):
            mw_arr = np.append(mw_arr, np.array(self.retrive_param_from_DB(cas, 'MW'), dtype=np.float32))
            a_num = np.append(a_num, np.array(self.retrive_param_from_DB(cas, 'ANum'), dtype=np.float32))
            ac_ = self.retrive_param_from_DB(cas, 'AC')
            number_array = list(map(int, ac_))
            ac_arr = np.array(number_array)
            ac_sum = np.sum(ac_arr * self.param['atomic constants'])
            ac_sum_arr = np.append(ac_sum_arr, ac_sum)
            vdw_arr = np.append(vdw_arr, np.array(self.retrive_param_from_DB(cas, 'DVAMVDW'), dtype=np.float32))

        A = np.sum(ac_sum_arr * co_molefrac) / np.sum(vdw_arr * co_molefrac * 0.6022)

        tc = A * self.compute_solid_heat_capacity_Bicerano(param, T,
                                                           co_molefrac) * self.compute_solid_density_Askadskii_Matveev(
            param, T, co_molefrac) ** (4 / 3) / (np.sum(mw_arr * co_molefrac) / np.sum(a_num * co_molefrac)) ** (1 / 3)

        return tc / 100

    def compute_solid_Enthalpy_Bicerano(self, T, param):

        pass

    def compute_thermal_conductivity_TRAPP(self, T, mass_frac):

        return 0.148752177757864 + 0.000983609776862362 * (
                T - 273.15) - 0.0407279003983628 * mass_frac + 0.000282659734914412 * (T - 273.15) * mass_frac
