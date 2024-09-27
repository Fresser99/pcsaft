import pandas as pd
import numpy as np
from param import *
from pcsaft import *
import math

df = pd.DataFrame([[10500., 10500, 10500, 1400, 1300, 1380, 0., 0., 0., 3.7, 2.8, 3.1, -92.1, -93.8, -92.3, -94.1,
                    -91.7, -94.2, 113.1, 121.1, 131.1, 121.1, 111.1, 108.1, 37.1, 2.9, 74.1, 73.2, 73.5]])


def calculate_mole_fraction(w_a):
    m_a = 56
    m_b = 50.487
    w_b = 1 - w_a
    x = (w_a / m_a) / ((w_a / m_a) + (w_b / m_b))
    return x


w_IB_R130 = df.values[0, 9] / 100
w_IB_R140 = df.values[0, 10] / 100
w_IB_R150 = df.values[0, 11] / 100

x_a_130 = calculate_mole_fraction(w_IB_R130)
x_b_130 = 1 - x_a_130
mole_frac_R130 = np.array([x_a_130, x_b_130])

x_a_140 = calculate_mole_fraction(w_IB_R140)
x_b_140 = 1 - x_a_140
mole_frac_R140 = np.array([x_a_140, x_b_140])

x_a_150 = calculate_mole_fraction(w_IB_R150)
x_b_150 = 1 - x_a_150
mole_frac_R150 = np.array([x_a_150, x_b_150])

method = PCSAFT()
param = Param
param.m = np.array(
    [method.retrive_param_from_DB('115-11-7', 'PCFTM'), method.retrive_param_from_DB('74-87-3', 'PCFTM')],
    dtype=np.float32)
param.e = np.array(
    [method.retrive_param_from_DB('115-11-7', 'PCFTU'), method.retrive_param_from_DB('74-87-3', 'PCFTU')],
    dtype=np.float32)
param.s = np.array(
    [method.retrive_param_from_DB('115-11-7', 'PCFTV'), method.retrive_param_from_DB('74-87-3', 'PCFTV')],
    dtype=np.float32)
param.k_ij = np.zeros([2, 2])
param.type = np.array([0, 0])
param.CPIG = np.array([[35463.9, 229574.8, 1166.14, 116871.6, 394.9803, 25, 1226.85],
                       [36220, 69810, 1805, 44470, 844.27, 298.15, 1500]])
param.HIGTREF = np.array([-17.1, -85.7])

param.CPDIP=np.array([[87680,217.1,-0.9153,0.002266,0,132.81,343.15],[107900,-330.13,0.808,0,0,175.43,303.15]])

param_IIR = Param
param_IIR.CAS = np.array(['seg-IB-R', 'seg-IP-R'])
param_IIR.type = np.array([1, 1])

T_130_top = df.values[0, 12] + 273.15
T_130_bot = df.values[0, 13] + 273.15
p_130_top = df.values[0, 19] * 1000
p_130_bot = df.values[0, 18] * 1000

T_140_top = df.values[0, 14] + 273.15
T_140_bot = df.values[0, 15] + 273.15
p_140_top = df.values[0, 21] * 1000
p_140_bot = df.values[0, 20] * 1000

T_150_top = df.values[0, 16] + 273.15
T_150_bot = df.values[0, 17] + 273.15
p_150_top = df.values[0, 23] * 1000
p_150_bot = df.values[0, 22] * 1000


def logarithmic_mean(a, b):
    if a <= 0 or b <= 0:
        raise ValueError("Both numbers must be positive")

    if a == b:
        return a

    return (b - a) / (math.log(b) - math.log(a))


T_130_log_mean = logarithmic_mean(T_130_top, T_130_bot)
p_130_log_mean = logarithmic_mean(p_130_top, p_130_bot)

T_140_log_mean = logarithmic_mean(T_140_top, T_140_bot)
p_140_log_mean = logarithmic_mean(p_140_top, p_140_bot)

T_150_log_mean = logarithmic_mean(T_150_top, T_150_bot)
p_150_log_mean = logarithmic_mean(p_150_top, p_150_bot)

tot_IB_mass_130 = df.values[0, 0] * df.values[0, 24] / 100
tot_IB_mass_140 = df.values[0, 1] * df.values[0, 24] / 100
tot_IB_mass_150 = df.values[0, 2] * df.values[0, 24] / 100

conver_IB_130 = (df.values[0, 0] * (df.values[0, 24] / 100 - df.values[0, 9] / 100) - df.values[0, 9] / 100 * df.values[
    0, 3] - df.values[0, 9] / 100 * df.values[0, 6]) / (
                        1 - df.values[0, 9] / 100)
conver_IB_140 = (df.values[0, 1] * (df.values[0, 24] / 100 - df.values[0, 10] / 100) - df.values[0, 10] / 100 *
                 df.values[0, 4] - df.values[0, 10] / 100 * df.values[0, 7]) / (
                        1 - df.values[0, 10] / 100)
conver_IB_150 = (df.values[0, 2] * (df.values[0, 24] / 100 - df.values[0, 11] / 100) - df.values[0, 11] / 100 *
                 df.values[0, 5] - df.values[0, 11] / 100 * df.values[0, 8]) / (
                        1 - df.values[0, 11] / 100)

conver_IP_130 = df.values[0, 0] * (df.values[0, 24] / 100 * df.values[0, 25] / 100)
conver_IP_140 = df.values[0, 1] * (df.values[0, 24] / 100 * df.values[0, 25] / 100)
conver_IP_150 = df.values[0, 2] * (df.values[0, 24] / 100 * df.values[0, 25] / 100)

mol_a_130 = conver_IB_130 / 56 / (conver_IB_130 / 56 + conver_IP_130 / 68)
mol_b_130 = 1 - mol_a_130

mol_a_140 = conver_IB_140 / 56 / (conver_IB_140 / 56 + conver_IP_140 / 68)
mol_b_140 = 1 - mol_a_140

mol_a_150 = conver_IB_150 / 56 / (conver_IB_150 / 56 + conver_IP_150 / 68)
mol_b_150 = 1 - mol_a_150

co_molefrac_130 = np.array([mol_a_130, mol_b_130])
co_molefrac_140 = np.array([mol_a_140, mol_b_140])
co_molefrac_150 = np.array([mol_a_150, mol_b_150])

density_liq_130 = method.compute_density(T_130_log_mean, p_130_log_mean, mole_frac_R130, 0, param) / 1000
density_liq_140 = method.compute_density(T_140_log_mean, p_140_log_mean, mole_frac_R140, 0, param) / 1000
density_liq_150 = method.compute_density(T_150_log_mean, p_150_log_mean, mole_frac_R150, 0, param) / 1000

# 聚合物密度
density_sol_130 = method.compute_solid_density_Askadskii_Matveev(param_IIR, T_130_log_mean, co_molefrac_130) * 1000
density_sol_140 = method.compute_solid_density_Askadskii_Matveev(param_IIR, T_140_log_mean, co_molefrac_140) * 1000
density_sol_150 = method.compute_solid_density_Askadskii_Matveev(param_IIR, T_150_log_mean, co_molefrac_150) * 1000

tot_mole_liq_130 = (tot_IB_mass_130 - conver_IB_130) / 56 + (
        df.values[0, 0] * (1 - df.values[0, 24] / 100) + df.values[0, 3] + df.values[0, 6]) / 50.487

tot_mole_liq_140 = (tot_IB_mass_140 - conver_IB_140) / 56 + (
        df.values[0, 1] * (1 - df.values[0, 24] / 100) + df.values[0, 4] + df.values[0, 7]) / 50.487

tot_mole_liq_150 = (tot_IB_mass_150 - conver_IB_150) / 56 + (
        df.values[0, 2] * (1 - df.values[0, 24] / 100) + df.values[0, 5] + df.values[0, 8]) / 50.487

tot_v_liq_130 = tot_mole_liq_130 / density_liq_130
tot_v_liq_140 = tot_mole_liq_140 / density_liq_140
tot_v_liq_150 = tot_mole_liq_150 / density_liq_150

tot_v_sol_130 = (conver_IB_130 + conver_IP_130) / (density_sol_130)
tot_v_sol_140 = (conver_IB_140 + conver_IP_140) / (density_sol_140)
tot_v_sol_150 = (conver_IB_150 + conver_IP_150) / (density_sol_150)

mass_frac_130 = np.array([1 - (conver_IB_130 + conver_IP_130) / (df.values[0, 0] + df.values[0, 3] + df.values[0, 6]),
                          (conver_IB_130 + conver_IP_130) / (df.values[0, 0] + df.values[0, 3] + df.values[0, 6])])
mass_frac_140 = np.array([1 - (conver_IB_140 + conver_IP_140) / (df.values[0, 1] + df.values[0, 4] + df.values[0, 7]),
                          (conver_IB_140 + conver_IP_140) / (df.values[0, 1] + df.values[0, 4] + df.values[0, 7])])
mass_frac_150 = np.array([1 - (conver_IB_150 + conver_IP_150) / (df.values[0, 2] + df.values[0, 5] + df.values[0, 8]),
                          (conver_IB_150 + conver_IP_150) / (df.values[0, 2] + df.values[0, 5] + df.values[0, 8])])
# 淤浆密度
density_slurry_130 = (df.values[0, 0] + df.values[0, 3] + df.values[0, 6]) / (tot_v_sol_130 + tot_v_liq_130)
density_slurry_140 = (df.values[0, 1] + df.values[0, 4] + df.values[0, 7]) / (tot_v_sol_140 + tot_v_liq_140)
density_slurry_150 = (df.values[0, 2] + df.values[0, 5] + df.values[0, 8]) / (tot_v_sol_150 + tot_v_liq_150)

cp_liq_130 = np.sum(
    method.compute_cpig(T_130_log_mean, param, np.array([0, 0]))*mole_frac_R130) / np.sum(np.array([56, 50.4875]) * np.array(
        [w_IB_R130, 1 - w_IB_R130])) / 1000
cp_liq_140 = np.sum(
    method.compute_cpig(T_140_log_mean, param, np.array([0, 0]))*mole_frac_R140) / np.sum(np.array([56, 50.4875]) * np.array(
        [w_IB_R140, 1 - w_IB_R140])) / 1000
cp_liq_150 = np.sum(
    method.compute_cpig(T_150_log_mean, param, np.array([0, 0]))*mole_frac_R150) / np.sum(np.array([56, 50.4875]) * np.array(
        [w_IB_R150, 1 - w_IB_R150])) / 1000

cp_liq_130_mix=np.sum(method.compute_liquid_heat_capacity(T_130_log_mean,param) * mole_frac_R130) / np.sum(np.array([56, 50.4875]) * np.array(
        [w_IB_R130, 1 - w_IB_R130])) / 1000
cp_liq_140_mix=np.sum(method.compute_liquid_heat_capacity(T_140_log_mean,param)*mole_frac_R140)/np.sum(np.array([56, 50.4875]) * np.array(
        [w_IB_R140, 1 - w_IB_R140])) / 1000
cp_liq_150_mix=np.sum(method.compute_liquid_heat_capacity(T_150_log_mean,param)*mole_frac_R150)/np.sum(np.array([56, 50.4875]) * np.array(
        [w_IB_R150, 1 - w_IB_R150])) / 1000

# 聚合物热容
cp_sol_130 = method.compute_solid_heat_capacity_Bicerano(param_IIR, T_130_log_mean, co_molefrac_130) / np.sum(
    np.array([56, 68]) * co_molefrac_130)
cp_sol_140 = method.compute_solid_heat_capacity_Bicerano(param_IIR, T_140_log_mean, co_molefrac_140) / np.sum(
    np.array([56, 68]) * co_molefrac_140)
cp_sol_150 = method.compute_solid_heat_capacity_Bicerano(param_IIR, T_150_log_mean, co_molefrac_150) / np.sum(
    np.array([56, 68]) * co_molefrac_150)

v_frac_phase_130 = np.array(
    [tot_v_liq_130 / (tot_v_sol_130 + tot_v_liq_130), tot_v_sol_130 / (tot_v_sol_130 + tot_v_liq_130)])
v_frac_phase_140 = np.array(
    [tot_v_liq_140 / (tot_v_sol_140 + tot_v_liq_140), tot_v_sol_140 / (tot_v_sol_140 + tot_v_liq_140)])
v_frac_phase_150 = np.array(
    [tot_v_liq_150 / (tot_v_sol_150 + tot_v_liq_150), tot_v_sol_150 / (tot_v_sol_150 + tot_v_liq_150)])

# 淤浆热容
cp_slurry_130 = np.sum(mass_frac_130 * np.array([cp_liq_130_mix, cp_sol_130]))
cp_slurry_140 = np.sum(mass_frac_140 * np.array([cp_liq_140_mix, cp_sol_140]))
cp_slurry_150 = np.sum(mass_frac_150 * np.array([cp_liq_150_mix, cp_sol_150]))

# hres_liq_130 = method.compute_hres(T_130_log_mean, density_liq_130, mole_frac_R130, param)
# hres_liq_140 = method.compute_hres(T_140_log_mean, density_liq_140, mole_frac_R140, param)
# hres_liq_150 = method.compute_hres(T_150_log_mean, density_liq_150, mole_frac_R150, param)
#
# hig_liq_arr_130 = method.compute_hig(param, T_130_log_mean)
# hig_liq_arr_140 = method.compute_hig(param, T_140_log_mean)
# hig_liq_arr_150 = method.compute_hig(param, T_150_log_mean)
#
# hig_liq_mix_130 = method.compute_hig_mix(hig_liq_arr_130, mole_frac_R130)
# hig_liq_mix_140 = method.compute_hig_mix(hig_liq_arr_140, mole_frac_R140)
# hig_liq_mix_150 = method.compute_hig_mix(hig_liq_arr_150, mole_frac_R150)

# 淤浆焓
# h_liq_130 = method.compute_Enthalpy(hig_liq_mix_130, hres_liq_130)
# h_liq_140 = method.compute_Enthalpy(hig_liq_mix_140, hres_liq_140)
# h_liq_150 = method.compute_Enthalpy(hig_liq_mix_150, hres_liq_150)

param3 = Param
param3.m = np.array([0.0238 * 10000])
param3.e = np.array([223])
param3.s = np.array([4.088])
param3.type = np.array([1])
param3.k_ij = np.array([0.])
param3.CPIG = np.array([-66039, 715.84, -0.7804, 0.0003255, 0, 0, 280, 1000, 36029.2, 0.142427, 2.244683])
param3.HIGTREF = np.array([0.])

# den_sol_130 = method.compute_density(T_130_log_mean, p_130_log_mean, np.array([1.]), 0, param3)
# hres_sol_130 = method.compute_hres(T_130_log_mean, den_sol_130, np.array([1.]), param3)
# hres_sol_cor_130 = method.compute_hres_correction(hres_sol_130, 10000 / 56)
# hig_sol_130 = method.compute_hig(param3, T_130_log_mean)
# # 聚合物焓
# hig_sol_mix_130 = method.compute_hig_mix(hig_sol_130, np.array([1.]))
# h_sol_130 = hres_sol_cor_130 + hig_sol_mix_130

# den_sol_140 = method.compute_density(T_140_log_mean, p_140_log_mean, np.array([1.]), 0, param3)
# hres_sol_140 = method.compute_hres(T_140_log_mean, den_sol_140, np.array([1.]), param3)
# hres_sol_cor_140 = method.compute_hres_correction(hres_sol_140, 10000 / 56)
# hig_sol_140 = method.compute_hig(param3, T_140_log_mean)
# # 聚合物焓
# hig_sol_mix_140 = method.compute_hig_mix(hig_sol_140, np.array([1.]))
# h_sol_140 = hres_sol_cor_140 + hig_sol_mix_140
#
# den_sol_150 = method.compute_density(T_150_log_mean, p_150_log_mean, np.array([1.]), 0, param3)
# hres_sol_150 = method.compute_hres(T_150_log_mean, den_sol_150, np.array([1.]), param3)
# hres_sol_cor_150 = method.compute_hres_correction(hres_sol_150, 10000 / 56)
# hig_sol_150 = method.compute_hig(param3, T_150_log_mean)
# # 聚合物焓
# hig_sol_mix_150 = method.compute_hig_mix(hig_sol_150, np.array([1.]))
# h_sol_150 = hres_sol_cor_150 + hig_sol_mix_150

# 聚合物热导率
tc_sol_130 = method.compute_solid_thermal_conductivity_Askadskii_Matveev(param_IIR, co_molefrac_130, T_130_log_mean)
tc_sol_140 = method.compute_solid_thermal_conductivity_Askadskii_Matveev(param_IIR, co_molefrac_140, T_140_log_mean)
tc_sol_150 = method.compute_solid_thermal_conductivity_Askadskii_Matveev(param_IIR, co_molefrac_150, T_150_log_mean)

tc_liq_130 = method.compute_thermal_conductivity_TRAPP(T_130_log_mean, df.values[0, 9] / 100)
tc_liq_140 = method.compute_thermal_conductivity_TRAPP(T_140_log_mean, df.values[0, 10] / 100)
tc_liq_150 = method.compute_thermal_conductivity_TRAPP(T_150_log_mean, df.values[0, 11] / 100)

tc_slurry_130 = np.sum(np.array([tc_liq_130, tc_sol_130]) * v_frac_phase_130)
tc_slurry_140 = np.sum(np.array([tc_liq_140, tc_sol_140]) * v_frac_phase_140)
tc_slurry_150 = np.sum(np.array([tc_liq_150, tc_sol_150]) * v_frac_phase_150)

#
T_135 = df.values[0, 26] + 273.15
T_145 = df.values[0, 27] + 273.15
T_155 = df.values[0, 28] + 273.15

density_sol_135 = method.compute_solid_density_Askadskii_Matveev(param_IIR, T_135, co_molefrac_130) * 1000
density_sol_145 = method.compute_solid_density_Askadskii_Matveev(param_IIR, T_145, co_molefrac_140) * 1000
density_sol_155 = method.compute_solid_density_Askadskii_Matveev(param_IIR, T_155, co_molefrac_150) * 1000

cp_sol_135 = method.compute_solid_heat_capacity_Bicerano(param_IIR, T_135, co_molefrac_130) / np.sum(
    np.array([56, 68] * co_molefrac_130))
cp_sol_145 = method.compute_solid_heat_capacity_Bicerano(param_IIR, T_145, co_molefrac_140) / np.sum(
    np.array([56, 68] * co_molefrac_140))
cp_sol_155 = method.compute_solid_heat_capacity_Bicerano(param_IIR, T_155, co_molefrac_150) / np.sum(
    np.array([56, 68] * co_molefrac_150))

tc_sol_135 = method.compute_solid_thermal_conductivity_Askadskii_Matveev(param_IIR, co_molefrac_130, T_135)
tc_sol_145 = method.compute_solid_thermal_conductivity_Askadskii_Matveev(param_IIR, co_molefrac_140, T_145)
tc_sol_155 = method.compute_solid_thermal_conductivity_Askadskii_Matveev(param_IIR, co_molefrac_150, T_155)

Tg_130 = method.compute_glassify_temperature_Askadskii_Matveev(param_IIR, co_molefrac_130)
Tg_140 = method.compute_glassify_temperature_Askadskii_Matveev(param_IIR, co_molefrac_140)
Tg_150 = method.compute_glassify_temperature_Askadskii_Matveev(param_IIR, co_molefrac_150)
