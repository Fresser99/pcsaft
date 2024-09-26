import numpy as np
import math
from pcsaft import *
from param import *
import redis
import requests


def main():
    # 读取数据
    datasource_config = {
        'host': '10.0.1.211',
        'port': '30379',
        'db': '1',
        'password': '',
        'start': -1,
        'end': -1,
        'group_id': 1824350503343808512,
    }

    def read_redis(host: str, port: str, db: str, group_id: str, tag_list: list, password: str = '', start: int = 0,
                   end: int = -1):
        """从 redis 读取指定位号数据
        :param host: redis 地址
        :param port: 端口，通常为 30379
        :param db:分区
        :param password: 密码，默认为空
        :param group_id: 位号组 id
        :param tag_list: 位号列表
        :param start: 列表起始位置
        :param end:: 列表结束位置
        :return:
        """
        pool = redis.ConnectionPool(host=host, port=port, db=db, password=password, decode_responses=True)
        client = redis.Redis(connection_pool=pool)
        ls = [client.lrange(f'data:tag_value:group:{group_id}:{tag_id}', start, end) for tag_id in
              tag_list]
        df = pd.concat([pd.DataFrame([eval(eval(i)) for i in tag]) for tag in ls])
        if len(df) == 0:
            assert False, (', '.join(tag_list) + '不存在')
        del df['tag_time']
        df.rename(columns={'app_time': 'tag_time'}, inplace=True)
        df = df[['tag_time', 'tag_name', 'tag_value']]
        df['tag_value'] = pd.to_numeric(df['tag_value'], errors='coerce')
        df.drop_duplicates(subset=['tag_time', 'tag_name'], keep='first', inplace=True)  # 去除重复
        df = pd.pivot(df, values='tag_value', index='tag_time', columns='tag_name')
        df.dropna(axis=0, how='any', inplace=True)  # 去掉空值
        df.index = pd.to_datetime(df.index)
        df = df.sort_values('tag_time')
        client.close()
        return df

    tag_list = ['FIC_A13002', 'FIC_A14002', 'FIC_A15002', 'FIC_A13001A', 'FIC_A14001A',
                'FIC_A15001A', 'FIC_A13001B', 'FIC_A14001B', 'FIC_A15001B', 'AI_A13502', 'AI_A14502', 'AI_A15502',
                'TIC_A13004',
                'TI_A13001', 'TIC_A14004', 'TI_A14001', 'TIC_A15004', 'TI_A15001', 'PI_A13001', 'PI_A13008',
                'PI_A14001',
                'PI_A14008', 'PI_A15001', 'PI_A15008', '100_Oper.K8', '100_Oper.K9']

    df = read_redis(**datasource_config, tag_list=tag_list)

    def calculate_mole_fraction(w_a):
        m_a = 56
        m_b = 50.487
        w_b = 1 - w_a
        x = (w_a / m_a) / ((w_a / m_a) + (w_b / m_b))
        return x

    w_IB_R130 = df.values[9] / 100
    w_IB_R140 = df.values[10] / 100
    w_IB_R150 = df.values[11] / 100

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
    param.type = np.array([0])
    param.CPIG = np.array([[35463.9, 229574.8, 1166.14, 116871.6, 394.9803, 25, 1226.85],
                           [36220, 69810, 1805, 44470, 844.27, 298.15, 1500, 2.4161]])
    param.HIGTREF = np.array([-17.1, -85.7])

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

    conver_IB_130 = (df.values[0, 0] * (1 - df.values[0, 25]) + df.values[0, 3] + df.values[0, 6]) * df.values[0, 9] / (
            1 - df.values[0, 9])
    conver_IB_140 = (df.values[0, 1] * (1 - df.values[0, 25]) + df.values[0, 4] + df.values[0, 7]) * df.values[
        0, 10] / (
                            1 - df.values[0, 10])
    conver_IB_150 = (df.values[0, 2] * (1 - df.values[0, 25]) + df.values[0, 5] + df.values[0, 8]) * df.values[
        0, 11] / (
                            1 - df.values[0, 11])

    conver_IP_130=df.values[0,0]*

    mol_a_130 = df.values[0, 25] / 100 * 56 / (df.values[0, 25] * 56 + 68.1185)
    mol_b_130 = 1 - mol_a_130

    co_molefrac_130 = np.array([mol_b_130, mol_a_130])
    density_liq_130 = method.compute_density(T_130_log_mean, p_130_log_mean, mole_frac_R130, 0, param)
    # 聚合物密度
    density_sol_130 = method.compute_solid_density_Askadskii_Matveev(param_IIR, T_130_log_mean, co_molefrac_130)
    tot_mole_liq_130 = (tot_IB_mass_130 - conver_IB_130) / 56 + (
            df.values[0, 0] * (1 - df.values[0, 24]) + df.values[0, 3] + df.values[0, 6]) / 50.487
    tot_v_liq_130 = tot_mole_liq_130 / density_liq_130
    tot_v_sol_130 = conver_IB_130 / (density_sol_130 * 1000)
    # 淤浆密度
    density_slurry = (df.values[0, 0] + df.values[0, 3] + df.values[0, 6]) / (tot_v_sol + tot_v_liq)
    cp_liq = np.sum(method.compute_cpig(T_130_log_mean, param, np.array([0, 0])) * mole_frac)
    # 聚合物热容
    cp_sol = method.compute_solid_heat_capacity_Bicerano(param_IIR, T_130_log_mean, co_molefrac)

    v_frac_phase = np.array([tot_v_liq / (tot_v_sol + tot_v_liq), tot_v_sol / (tot_v_sol + tot_v_liq)])
    # 淤浆热容
    cp_slurry = np.sum(v_frac_phase * np.array([cp_liq, cp_sol]))

    hres_liq = method.compute_hres(T_130_log_mean, density_liq, mole_frac, param)
    hig_liq_arr = method.compute_hig(param, T_130_log_mean)
    hig_liq_mix = method.compute_hig_mix(hig_liq_arr, mole_frac)
    # 淤浆焓
    h_liq = method.compute_Enthalpy(hig_liq_mix, hres_liq)

    param3 = Param
    param3.m = np.array([0.0238 * 10000])
    param3.e = np.array([223])
    param3.s = np.array([4.088])
    param3.type = np.array([1])
    param3.k_ij = np.array([0.])
    param3.CPIG = np.array([-66039, 715.84, -0.7804, 0.0003255, 0, 0, 280, 1000, 36029.2, 0.142427, 2.244683])
    param3.HIGTREF = np.array([0.])

    den_sol = method.compute_density(T_130_log_mean, p_130_log_mean, np.array([1.]), 0, param3)
    hres_sol = method.compute_hres(T_130_log_mean, den_sol, np.array([1.]), param3)
    hres_sol_cor = method.compute_hres_correction(hres_liq, 10000 / 56)
    hig_sol = method.compute_hig(param3, T_130_log_mean)
    # 聚合物焓
    hig_sol_mix = method.compute_hig_mix(hig_sol, np.array([1.]))
    h_sol = hres_sol + hig_sol_mix
    # 聚合物热导率
    tc_sol = method.compute_solid_thermal_conductivity_Askadskii_Matveev(param_IIR, co_molefrac, T_130_log_mean)

    # 将结果保存到自定义位号中
    url = "http://10.0.1.211:31013/api/tag-value/writeTagValues"

    headers = {
        'Connection': 'keep-alive',
        'Content-Type': 'application/json', }
    request_info = {
        "data": {
            "values": {
                "R130 slurry density": density_slurry,
                "R130 Polymer density": density_sol,
                "R130 slurry heat capacity": cp_slurry,
                "R130 polymer heat capacity": cp_sol,
                "R130 slurry enthalpy": h_liq,
                "R130 polymer enthalpy": h_sol,
                "R130 polymer thermal conductivity": tc_sol
            }
        }, "requestBase": {
            "page": "1-10", "sort": "-createTime"}
    }
    response = requests.post(url, json=request_info, headers=headers)


if __name__ == '__main__':
    main()
