from functools import lru_cache
from typing import Tuple
from scipy.optimize import brentq, differential_evolution
import numpy as np


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
            'conv': 7242.702976750923
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

            num_pts = 50  # 增加点数以提高精度

            rho_guesses = np.logspace(-13, np.log10(0.7405), num_pts)
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

    @lru_cache(maxsize=128)
    def _solve_density(self, rho_low: float, rho_up: float, T: float, p: float, comp: np.ndarray,
                       param) -> float:
        """使用Brent方法求解密度"""
        rho_low = self.compute_molar_density(rho_low, T, len(comp), comp, param)
        rho_up = self.compute_molar_density(rho_up, T, len(comp), comp, param)

        def f(rho):
            return self.solve_rho(rho, T, p, comp, param)

        return brentq(f, rho_low, rho_up, xtol=1e-8, maxiter=200)

    def _solve_multiple_roots(self, x_low, x_up, T:float, p: float, comp: np.ndarray,
                              param) -> float:
        """使用差分进化算法处理多个根的情况"""
        bounds = list(zip(x_low, x_up))
        result = differential_evolution(lambda x: self.compute_pcsaft_gibbs(T, x[0], comp, param),
                                        bounds, popsize=20, tol=1e-8)
        return result.x[0]

    def _global_optimization(self, T: float, p: float, comp: List[float], param: Dict) -> float:
        """当没有找到根时使用全局优化"""
        bounds = [(1e-8, 7.4)]
        result = differential_evolution(
            lambda x: abs(self.compute_pressure_residual(x[0], T, p, comp, param)),
            bounds, popsize=20, tol=1e-8)
        return self.compute_pcsaft_molar_density(result.x[0], T, len(comp), comp, param)

    def _find_root_intervals(self, rho_guesses: np.ndarray, T: float, p: float, comp: List[float], param: Dict) -> \
    Tuple[List[float], List[float]]:
        """找出可能的根区间"""
        x_low, x_up = [], []
        P_err_prev = self.compute_pressure_residual(rho_guesses[0], T, p, comp, param)

        for rho_guess in rho_guesses[1:]:
            P_err = self.compute_pressure_residual(rho_guess, T, p, comp, param)

            if P_err_prev * P_err < 0:
                x_low.append(rho_guesses[rho_guesses < rho_guess][-1])
                x_up.append(rho_guess)

            P_err_prev = P_err

        return x_low, x_up

    def _validate_inputs(self, T: float, p: float, comp: List[float], phase: int, param: Dict):
        """验证输入参数"""
        if T <= 0 or p <= 0:
            raise ValueError("温度和压力必须为正值")
        if not (0 <= phase <= 1):
            raise ValueError("相态必须为0（液相）或1（气相）")
        if not (0.99 <= sum(comp) <= 1.01):
            raise ValueError("组分之和必须接近1")
