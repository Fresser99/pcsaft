
# PCSAFT 使用说明

作者：ZFY
日期：9月22日

## 部分测试结果
### 低聚物摩尔密度（对比Aspen）
<img width="415" alt="image" src="https://github.com/user-attachments/assets/9cd74dd6-53bb-46e1-b559-c2f61df63ae2" />

### 混合物系

#### 流股条件

| 组分 | 含量 |  
|------|------|  
| TICL4 | 0.036902158 |  
| TEA | 0.091969849 |  
| C2H4 | 76.98073312 |  
| H2 | 1.082327566 |  
| HEXANE | 719.2184100 |  
| HDPE | 0.150378939 |  
| 温度 | 160 K |  
| 压力 | 200 atm |  
| DPN | 993.43 |  

| 组分 | 含量 |  
|------|------|  
| TICL4 | 3.76E-05 |  
| TEA | 9.38E-05 |  
| C2H4 | 0.104064454 |  
| H2 | 0.000164843 |  
| HEXANE | 0.74526035 |  
| HDPE | 0.150378939 |  
| 温度 | 160 K |  
| 压力 | 200 atm |  
| DPN | 1174.19 |  


| 物性 | sup-CSTR | Aspen | 相对偏差 [%] |  
|------|----------|-------|-------------|  
| 摩尔密度 | 7908.830319 | 7908.96487 | 0.001701277 |  
| 压缩因子 | 0.711496907 | 0.711477852 | 0.002678255 |  
| 焓偏离 | -20147.78624 | -20146.95233 | -0.004138999 |  
| 熵偏离 | -66.89850028 | -66.89569679 | -0.004190655 |  
| 吉布斯偏离 | -942806.7513 | -942806.9152 | -1.73808E-05 |  
| **逸度系数** |  |  |  |  
| TICL4 | 2.01E-14 | 2.01E-14 | 0.000374278 |  
| TEA | 2.01E-14 | 2.01E-14 | 0.000374278 |  
| C2H4 | 1.112936799 | 1.112937281 | 4.32878E-05 |  
| H2 | 6.508317556 | 6.50832283 | 8.10328E-05 |  
| HEXANE | 0.102929471 | 0.102929477 | 6.27465E-06 |  
| HDPE | 0 | 2.67E-109 | 0 |  


| 物性 | sup-CSTR | Aspen | 相对偏差 [%] |  
|------|----------|-------|-------------|  
| 摩尔密度 | 8255.046764 | 8255.205643 | 0.00192463 |  
| 压缩因子 | 0.681683306 | 0.681636967 | 0.006797834 |  
| 焓偏离 | -19936.81131 | -19935.18659 | -0.008149365 |  
| 熵偏离 | -63.84696831 | -63.84166453 | -0.008307017 |  
| **逸度系数** |  |  |  |  
| TICL4 | 2.14E-14 | 2.14E-14 | 0.000466137 |  
| TEA | 2.14E-14 | 2.14E-14 | 0.000466137 |  
| C2H4 | 1.201831619 | 1.201832592 | 8.09234E-05 |  
| H2 | 7.138994454 | 7.139006704 | 0.000171597 |  
| HEXANE | 0.110429903 | 0.110429919 | 1.45005E-05 |  
| HDPE | 0 | 2.67E-109 | 0 |

### 闪蒸计算结果

#### 流股条件
| 组分 | 含量 |  
|------|------|  
| TICL4 | 3.76E-05 |  
| TEA | 9.38E-05 |  
| C2H4 | 0.104064454 |  
| H2 | 0.000164843 |  
| HEXANE | 0.74526035 |  
| HDPE | 0.150378939 |  
| 温度 | 160 K |  
| 压力 | 10 atm |  
| DPN | 1174.19 |  

#### 对比结果（气相摩尔组成和液相摩尔组成 VS Aspen)

| 组分 | ASPEN 液相摩尔组成 | Sup-CSTR液相摩尔组成 | 相对偏差 [%] |  
|------|-------------------|---------------------|-------------|  
| TICL4 | 2.25E-04 | 0.000225272 | 0.00040602 |  
| TEA | 5.61E-04 | 0.000561437 | 0.000407497 |  
| C2H4 | 0.001268745 | 0.001268874 | 0.010169121 |  
| H2 | 2.84E-07 | 2.84E-07 | 0.010709016 |  
| HEXANE | 0.097697726 | 0.097709037 | 0.011576711 |  
| HDPE | 9.00E-01 | 0.900235095 | 0.001270485 |  

| 组分 | ASPEN 汽相摩尔组成 | Sup_CSTR汽相摩尔组成 | 相对偏差 [%] |  
|------|-------------------|----------------------|-------------|  
| TICL4 | 2.24E-15 | 2.24E-15 | 0.05129298 |  
| TEA | 5.59E-15 | 5.58E-15 | 0.051294458 |  
| C2H4 | 0.124679173 | 0.124679442 | 0.000215933 |  
| H2 | 0.000197844 | 0.000197844 | 0.000234491 |  
| HEXANE | 0.875122983 | 0.875122713 | 3.08172E-05 |  
| HDPE | 2.40E-109 | 0 | 0 |  

#### 气液平衡常数（VS Aspen）

| 组分 | ASPEN 气液平衡常数 | sup-CSTR气液平衡常数 | 相对偏差 |  
|------|---------------------|----------------------|---------|  
| TICL4 | 9.95E-12 | 9.94E-12 | 0.050816362 |  
| TEA | 9.95E-12 | 9.94E-12 | 0.050816362 |  
| C2H4 | 98.2696581 | 98.25989071 | 0.009940359 |  
| H2 | 695.6085917 | 695.5358307 | 0.010461145 |  
| HEXANE | 8.957453562 | 8.956415268 | 0.011592743 |  
| HDPE | 0 | 0 | 0 |

## 算法说明

### 函数汇总表

| 函数名称                                                         | 描述                         |
|------------------------------------------------------------------|----------------------------|
| `__init__(self)`                                                 | 初始化 PCSAFT 类的实例。           |
| `param_init(self)`                                               | 设置 PC-SAFT 模型的默认参数。        |
| `compute_z(self, rho, T, comp, param)`                           | 计算压缩因子 `Z`。                |
| `_compute_hs_terms(self, d, zeta)`                               | 计算硬球径向分布函数 `g_hs` 及其密度导数。  |
| `_compute_zhs(self, zeta)`                                       | 计算硬球压缩因子 `Zhs`。            |
| `_compute_ab(self, m_avg)`                                       | 计算色散相互作用的系数 `a` 和 `b`。     |
| `_compute_I_terms(self, eta3, a, b)`                             | 计算色散能量的积分系数。               |
| `_compute_C_terms(self, m_avg, eta3)`                            | 计算色散相关中的附加修正项 `C1` 和 `C2`。 |
| `compute_hres(self, T, rho, comp, param)`                        | 计算焓偏离 `hres`。              |
| `compute_sres(self, T, rho, comp, param)`                        | 计算熵偏离 `Sres`。              |
| `compute_Ares(self, T, rho, comp, param)`                        | 计算亥姆霍兹能偏离 `Ares`。          |
| `compute_gres(self, T, rho, comp, param)`                        | 计算吉布斯自由能偏离 `Gres`。         |
| `compute_density(self, T, p, comp, phase, param)`                | 计算给定温度和压力下的密度 `rho`。       |
| `_solve_density(self, rho_low, rho_up, T, p, comp, param)`       | 使用 Brent 方法求解密度。           |
| `_solve_multiple_roots(self, x_low, x_up, T, p, comp, param)`    | 处理具有多个密度根的场景。              |
| `_global_optimization(self, T, p, comp, param)`                  | 应用全局优化策略以找到真实物理状态的根。       |
| `compute_molar_density(self, rho_low, T, n, comp, param)`        | 将初始密度猜测转换为摩尔密度。            |
| `compute_pressure_residual(self, rho_guess, T, p, comp, param)`  | 计算密度猜测的压力残差。               |
| `compute_den_correction(self, rho, polymer_zeroth_mole_frac, mn)` | 计算并返回修正后的密度。               |
| `compute_z_correction(self, z, polymer_first_mole_frac, polymer_zeroth_mole_frac)` | 计算并返回修正后的压缩因子。             |
| `compute_hres_correction(self, hres, mn)`                        | 计算并返回修正后的焓偏离。              |
| `compute_hig(self, param, T)`                                    | 计算给定温度的理想气体焓。              |
| `compute_hig_mix(self, hig_arr, comp)`                           | 计算混合物的理想气体焓。               |
| `compute_cpig(self, T, param, flag)`                             | 计算混合物的理想气体定压比热容。           |
| `compute_Enthalpy(self, hig, hres)`                              | 计算混合物的焓。                   |
| `retrive_param_from_DB(self, CAS, param_name)`                   | 从数据库检索特定参数。                |
| `compute_glassify_temperature_Askadskii_Matveev(self, param, co_molefrac)` | 计算共聚体系的玻璃化温度。              |
| `compute_solid_density_Askadskii_Matveev(self, param, T, co_molefrac)` | 计算固相聚合物密度。                 |
| `compute_solid_heat_capacity_Bicerano(self, param, T, co_molefrac)` | 计算固相聚合物的热容。                |
| `compute_solid_thermal_conductivity_Askadskii_Matveev(self, param, co_molefrac, T)` | 计算固相聚合物的热导率。               |

## 使用方法
### 安装

使用git命令下载至本地文件夹
``` bash
git clone https://github.com/Fresser99/pcsaft.git
```

### 使用

#### 普通体系
``` python
from typing import NamedTuple
import numpy as np
from pcsaft import *
from param import *

#实例化PCSAFT类
method = PCSAFT()
#实例化参数类
param = Param()
param.m = np.array([1, 2, 3, 4])
param.e = np.array([200, 300, 100, 180])
param.s = np.array([2, 3, 4, 5])
param.k_ij = np.zeros([4, 4])
comp = np.array([0.3, 0.2, 0.1, 0.4])
T = 298.15 #Kelvin
p = 101325 #Pa
#计算摩尔密度
rhotest=method.compute_density(T,p,comp,0,param)
print(rhotest)
#计算压缩因子
ztest=method.compute_z(rhotest,T,comp,param)
print(ztest)
```

#### 聚合物体系
``` python
from typing import NamedTuple
import numpy as np
from pcsaft import *
from param import *

#实例化PCSAFT类
method = PCSAFT()
#实例化参数类

#注意：这里假定第四个组分为聚合物，则其特征链尺寸参数(m)则应由 r*Mn计算
#其中 r 为 特征链尺寸比，Mn 为该聚合物的数均分子量: Mn=dpn*MW, dpn为数均聚合度
#MW为构成聚合物单体的分子量
param = Param()

#假定聚合物为聚乙烯
dpn=1174.8
MW=28.0538
r=0.04132

param.m = np.array([1, 2, 3, r*dpn*MW])

#其余参数不变
param.e = np.array([200, 300, 100, 180])
param.s = np.array([2, 3, 4, 5])
param.k_ij = np.zeros([4, 4])
comp = np.array([0.3, 0.2, 0.1, 0.4])

#定义以聚合物0阶矩为基准的摩尔分率
comp_0=np.array([0.,0.,0.,0.])

T = 298.15 #Kelvin
p = 101325 #Pa
#计算摩尔密度
rhotest=method.compute_density(T,p,comp,0,param)
#将直接计算得到的摩尔密度和0阶矩摩尔分率以及聚合物的数均聚合度传入
#compute_den_correction
rhomix=method.compute_den_correction(rhotest,comp_0,dpn)
print(rhomix)
```

## 函数说明

### `compute_z(self, rho: float, T: float, comp: np.ndarray, param)`

- **参数：**
  - `rho`：混合物的密度（**注意：该密度因为由`compute_density`直接计算得到的密度值**）
  - `T`：温度（开尔文）。
  - `comp`：混合物中组分的摩尔分数（**注意：如果体系中包含聚合物组分时，使用的应该为以聚合物0阶矩摩尔量为基准计算的摩尔分率**）
  - `param`：特定于组分的参数，通常包含段数（`s`）、能量（`e`）和尺寸参数（`m`）或（`r`）。
  
- **返回：** 压缩因子 `Z`。
- **算法：**
  - 根据密度和温度计算段直径 `d` 和填充分数 `zeta`。
  - 使用硬球方程计算辅助项 `g_hs` 和 `delta_ghs_rho`。
  - 确定各种交互参数和平均值。
  - 计算硬链贡献（`Zhc`）和色散项（`Zdisp`）。
  - 返回压缩因子 `Z`，它是硬链、色散和理想项的总和。

### `_compute_hs_terms(self, d: np.ndarray, zeta: np.ndarray)`

- **参数：**
  - `d`：段直径数组。
  - `zeta`：填充分数数组。
  
- **返回：** `g_hs`（硬球径向分布函数）和 `delta_ghs_rho`（`g_hs` 的密度导数）的元组。
- **算法：**
  - 基于段直径和填充分数计算 `g_hs` 及其导数。

### `_compute_zhs(self, zeta: np.ndarray)`

- **参数：**
  - `zeta`：填充分数数组。
  
- **返回：** 硬球压缩因子 `Zhs`。
- **算法：**
  - 使用基于 `zeta` 的解析表达式计算 `Zhs`。

### `_compute_ab(self, m_avg: float)`

- **参数：**
  - `m_avg`：平均段数。
  
- **返回：** `a` 和 `b` 的元组，它们是色散相互作用的系数。
- **算法：**
  - 使用段数的平均值和预定义的系数数组计算系数 `a` 和 `b`。

### `_compute_I_terms(self, eta3: float, a: np.ndarray, b: np.ndarray)`

- **参数：**
  - `eta3`：三阶填充分数。
  - `a`：用于计算的系数数组。
  - `b`：用于计算的系数数组。
  
- **返回：** 积分项 `I2`、`deltaI1_det` 和 `deltaI2_det`。
- **算法：**
  - 使用 `eta3`、`a` 和 `b` 计算色散能量的积分系数。

### `_compute_C_terms(self, m_avg: float, eta3: float)`

- **参数：**
  - `m_avg`：平均段数。
  - `eta3`：三阶填充分数。
  
- **返回：** 元组 `C1` 和 `C2`，它们是色散相关中的附加修正项。
- **算法：**
  - 使用 `m_avg` 和填充分数计算修正项，以根据流体相和段相互作用调整色散贡献。

## 密度和压力计算

### `compute_density(self, T, p, comp, phase, param)`

- **参数：**
  - `T`：温度（开尔文）。
  - `p`：压力（帕斯卡）。
  - `comp`：摩尔分数。
  - `phase`：相态标志，0 表示液相，1 表示气相。
  - `param`：相互作用参数。
  
- **返回：** 对应于给定 `T` 和 `p` 的密度 `rho`。
- **算法：**
  - 验证输入并根据相态和可用根选择适当的方法来查找密度。
  - 根据确定的根区间，使用 `_solve_density`、`_solve_multiple_roots` 和 `_global_optimization` 等方法获得精确结果。

### `_solve_density(self, rho_low, rho_up, T, p, comp, param)`

- **参数：**
  - `rho_low`：密度下限。
  - `rho_up`：密度上限。
  - 其他参数如前所述。
  
- **返回：** 最优密度 `rho`。
- **算法：**
  - 在指定的区间内使用 Brent 方法求解密度，其中符号变化表示根。

### `_solve_multiple_roots(self, x_low, x_up, T, p, comp, param)`

- **参数：** 类似于 `_solve_density`，但用于多个根区间。
- **返回：** 密度 `rho`。
- **算法：**
  - 使用进化优化技术处理具有多个密度根的场景，以确定对应于最小吉布斯自由能的物理根。

### `_global_optimization(self, T, p, comp, param)`

- **参数：** 类似于之前的密度函数。
- **返回：** 全局优化的密度。
- **算法：**
  - 应用全局优化策略以找到可能是感兴趣的真实物理状态的根，当简单的区间方法失败时。

## 实用和支持功能
### `retrive_param_from_DB(self, CAS, param_name)`

- **参数：**
  - `CAS`：物质的 CAS 编号。
  - `param_name`：参数名称。
  
- **返回：** 特定参数的值。
- **算法：**
  - 从本地数据库文件 `db.csv` 中检索所需参数。

### `compute_molar_density(self, rho_low, T, n, comp, param)`

- **参数：**
  - `rho_low`：初始密度估计。
  
- **返回：** 摩尔密度。
- **算法：**
  - 使用组分属性将初始密度猜测转换为分子数密度。

### `compute_pressure_residual(self, rho_guess, T, p, comp, param)`

- **参数：** `rho_guess`、`T`、`p`、`comp`、`param` 如前所述。
- **返回：** 密度猜测的压力残差。
- **算法：**
  - 计算猜测的密度与预期压力的偏差。

## 焓、吉布斯自由能、熵和亥姆霍兹能偏离计算

### `compute_hres(self, T, rho, comp, param)`

- **参数：**
  - `T`：温度（开尔文）。
  - `rho`：密度。
  - `comp`：组分摩尔分数。
  - `param`：参数。
  
- **返回：** 焓偏离 `hres`。
- **算法：**
  - 计算压缩因子 `Z` 和温度导数 `dares_dt`。
  - 使用公式计算焓偏离，考虑温度、压缩因子和参数。

### `compute_sres(self, T, rho, comp, param)`

- **参数：**
  - `T`：温度（开尔文）。
  - `rho`：密度。
  - `comp`：组分摩尔分数。
  - `param`：参数。
  
- **返回：** 熵偏离 `Sres`。
- **算法：**
  - 计算亥姆霍兹能偏离 `ares` 和压缩因子 `Z`。
  - 使用公式计算熵偏离，考虑温度、亥姆霍兹能偏离、压缩因子和压力。

### `compute_Ares(self, T, rho, comp, param)`

- **参数：**
  - `T`：温度（开尔文）。
  - `rho`：密度。
  - `comp`：组分摩尔分数。
  - `param`：参数。
  
- **返回：** 亥姆霍兹能偏离 `Ares`。
- **算法：**
  - 计算并返回亥姆霍兹能偏离，考虑温度、密度和组分。

### `compute_gres(self, T, rho, comp, param)`

- **参数：**
  - `T`：温度（开尔文）
  - `rho`：密度。
  - `comp`：组分摩尔分数
  - `param`：参数
  
- **返回：** 吉布斯自由能偏离 `Gres`
- **算法：**
  - 计算亥姆霍兹能偏离 `Ares` 和熵偏离 `Sres`
  - 使用公式计算吉布斯自由能偏离，考虑温度、亥姆霍兹能偏离和熵偏离

## 修正计算方法

### `compute_den_correction(self, rho, polymer_zeroth_mole_frac, mn)`

- **参数：**
  - `rho`：密度
  - `polymer_zeroth_mole_frac`：聚合物的零阶摩尔分数
  - `mn`：聚合物的数均分子量
  
- **返回：** 修正后的密度
- **算法：**
  - 计算并返回修正后的密度，考虑聚合物的摩尔分数和平均分子量

### `compute_z_correction(self, z, polymer_first_mole_frac, polymer_zeroth_mole_frac)`

- **参数：**
  - `z`：压缩因子
  - `polymer_first_mole_frac`：聚合物的一阶摩尔分数
  - `polymer_zeroth_mole_frac`：聚合物的零阶摩尔分数
  
- **返回：** 修正后的压缩因子
- **算法：**
  - 计算并返回修正后的压缩因子，考虑聚合物的摩尔分数

### `compute_hres_correction(self, hres, mn)`

- **参数：**
  - `hres`：焓偏离。
  - `mn`：聚合物的数均分子量
  
- **返回：** 修正后的焓偏离
- **算法：**
  - 计算并返回修正后的焓偏离，考虑聚合物的平均分子量

### `compute_hig(self, param, T)`

- **参数：**
  - `param`：参数对象，包含物质的特性。
  - `T`：终止温度。
  
- **返回：** 理想气体焓 `hig`。
- **算法：**
  - 通过积分方法，计算从 298.15 K 到目标温度 `T` 的理想气体焓。

### `compute_hig_mix(self, hig_arr, comp)`

- **参数：**
  - `hig_arr`：理想气体焓数组。
  - `comp`：混合物组分摩尔分数。
  
- **返回：** 混合物的理想气体焓。
- **算法：**
  - 对每个组分的焓进行加权平均，获得混合物的理想气体焓。

### `compute_cpig(self, T, param, flag)`

- **参数：**
  - `T`：温度。
  - `param`：参数对象。
  - `flag`：标志数组，标识组分类型。
  
- **返回：** 理想气体定压比热容 `cpig`。
- **算法：**
  - 根据不同温度区间和组分类型，计算各个组分的定压比热容。

### `compute_Enthalpy(self, hig, hres)`

- **参数：**
  - `hig`：理想气体焓。
  - `hres`：焓偏离。
  
- **返回：** 混合物的焓。
- **算法：**
  - 合并理想气体焓和焓偏离，生成最终焓值。

## 固相聚合物物性计算
### `compute_glassify_temperature_Askadskii_Matveev(self, param, co_molefrac)`

- **参数：**
  - `param`：参数对象。
  - `co_molefrac`：共聚组成摩尔分率。
  
- **返回：** 玻璃化温度 `Tg`。
- **算法：**
  - 计算共聚体系的玻璃化温度，考虑摩尔分数和参数。

### `compute_solid_density_Askadskii_Matveev(self, param, T, co_molefrac)`

- **参数：**
  - `param`：参数对象。
  - `T`：温度。
  - `co_molefrac`：共聚组成摩尔分率
  
- **返回：** 固体密度。
- **算法：**
  - 计算固体密度，考虑温度、玻璃化温度和摩尔分数。

### `compute_solid_heat_capacity_Bicerano(self, param, T, co_molefrac)`

- **参数：**
  - `param`：参数对象。
  - `T`：温度。
  - `co_molefrac`：共聚组成摩尔分率。
  
- **返回：** 固体热容。
- **算法：**
  - 计算固体的热容，考虑温度和摩尔分数。

### `compute_solid_thermal_conductivity_Askadskii_Matveev(self, param, co_molefrac, T)`

- **参数：**
  - `param`：参数对象。
  - `co_molefrac`：共聚组成摩尔分率
  - `T`：温度。
  
- **返回：** 固体热导率。
- **算法：**
  - 计算固体的热导率，考虑温度、摩尔分数和参数。
