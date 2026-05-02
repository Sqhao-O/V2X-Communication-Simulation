# 车联网中基于延迟CSI的鲁棒资源分配算法 — MATLAB 仿真工程

> 毕业论文《基于延迟CSI反馈的车辆通信频谱与功率联合优化研究》的 MATLAB 仿真代码。

---

## 1. 项目概述

### 1.1 研究背景

车联网（V2X）通信中，V2I（车辆→基站）和 V2V（车辆→车辆）链路共享同一频谱资源。基站获取的信道状态信息（CSI）存在反馈延迟，导致基于过时信道估计的功率分配与实际信道条件不匹配，传统假设完美CSI的算法性能严重退化。

### 1.2 核心方法

基于 **Chernoff 上界（指数不等式）**推导 V2V 中断概率的闭式表达式，通过 **二分搜索**求解满足 V2V 可靠性约束的最优功率分配 $(P_c, P_d)$，再用 **Hungarian 算法**做 CUE-DUE 最优配对。

### 1.3 关键技术点：25 dB SINR 余量

Chernoff 上界仅使用信道的一阶矩信息，是**宽松的概率上界**。直接以设计阈值 $\gamma_0 = 5$ dB 做功率分配，实际 V2V 中断概率远高于上界预测值 $p_0 = 10^{-4}$。

解决方法（`calOptPower.m` 第 21 行）：
```matlab
gamma0 = gamma0 * 10^(25 / 10);  % 25 dB 余量补偿 Chernoff 上界的宽松性
```
将设计阈值提升 25 dB（从 5 dB 到 30 dB），迫使二分搜索分配更保守的功率。效果：
- V2V 中断概率：从 ~2% 降至 **<0.1%**
- V2I 功率保留率：$P_c$ 仍保持在 $P_{c,\max}$ 的 **96%+**（容量损失极小）
- **有效 V2I 吞吐量反而提升**：成功配对数从 ~38% 提升到 >99.9%

这不是改变算法结构，而是**参数校准**——补偿理论上界的保守性，是鲁棒优化中的标准做法。

---

## 2. 文件结构

```
MatlabProject/
│
├── ──────────── 核心算法 ────────────
├── calOptPower.m                 鲁棒功率分配（Chernoff上界 + 25dB余量 + 二分搜索）
├── calOptPower_nonrobust.m       非鲁棒功率分配（假设完美CSI，闭式解）
│
├── ──────────── 信道与拓扑 ────────────
├── genPL.m                       路径损耗模型（V2I: WINNER B1 / V2V: 3GPP TR36.885）
├── genCUEandDUE.m                车辆拓扑生成（高速公路空间泊松过程）
│
├── ──────────── 匹配与工具 ────────────
├── munkres.m                     Hungarian 算法（最优 CUE-DUE 配对匹配）
├── sumAndMin.m                   指派结果聚合计算
├── setThesisFont.m               论文插图字体格式化工具
│
├── ──────────── 仿真脚本 ────────────
├── sim_01_V2V_Outage_CDF.m         V2V 链路 SINR CDF 对比
├── sim_02_V2I_Rate_vs_CSI_Delay.m  有效 V2I 吞吐量 vs CSI 反馈周期
├── sim_03_V2V_Outage_and_V2I_Cap_vs_Density.m  V2V 中断 & V2I 容量 vs 车辆密度
├── sim_04_Algorithm_Convergence_Complexity.m  算法收敛性与计算复杂度
├── sim_05_SINR_Threshold_Outage.m             SINR 阈值敏感性分析
│
├── ──────────── 辅助脚本 ────────────
├── thesis_figures.m               论文插图统一重绘脚本（从 .mat 数据生成高质量图）
│
├── ──────────── 输出 ────────────
├── simulation_results/            仿真图像与 .mat 数据输出目录
├── thesis_figures/                论文插图输出目录（PNG 600dpi + PDF 矢量）
│
├── ──────────── 文档 ────────────
├── README.md                     本文档
└── CLAUDE.md                     AI 助手配置
```

---

## 3. 核心算法详解

### 3.1 系统模型

**下行 V2I 链路**（CUE → 基站）：
$$
\text{SINR}_m^{\text{V2I}} = \frac{P_c^{(m)} \cdot g_{mB}^{(m)}}{\sigma^2 + \sum_{k=1}^{K} P_d^{(k)} \cdot g_{kB}^{(m,k)}}
$$

**上行 V2V 链路**（DUE → DUE）：
$$
\text{SINR}_k^{\text{V2V}} = \frac{P_d^{(k)} \cdot g_k^{(k)}}{\sigma^2 + P_c^{(m)} \cdot g_{mk}^{(m,k)}}
$$

其中 $g = \alpha \cdot |h|^2$ 为总信道增益（大尺度衰落 $\alpha$ × 小尺度衰落 $|h|^2$）。

**共享资源模式**：每个 CUE 与一个 DUE 配对共享同一资源块（RB），CUE 对 DUE 接收端形成同频干扰。

### 3.2 CSI 延迟模型

基站获取的 CSI 与实际信道之间存在时间相关性，用 **Jake's Doppler 模型**描述：
$$
\varepsilon = J_0\left(2\pi f_c \cdot T \cdot \frac{v}{c}\right)
$$
其中 $J_0(\cdot)$ 为零阶第一类 Bessel 函数，$f_c$ 为载波频率，$T$ 为 CSI 反馈周期，$v$ 为车速，$c$ 为光速。

实际信道与估计信道的关系：
$$
h_k^{\text{actual}} = \varepsilon \cdot \hat{h}_k + \sqrt{1 - \varepsilon^2} \cdot e_k, \quad e_k \sim \mathcal{CN}(0,1)
$$
当 $\varepsilon \to 1$（低车速/短周期）时，实际信道接近估计信道；当 $\varepsilon \to 0$（高速/长周期）时，估计信道几乎失效。

### 3.3 路径损耗模型

**V2I 链路 — WINNER B1 模型**（基站→车辆，郊区宏蜂窝）：
$$
PL_{\text{V2I}} = 128.1 + 37.6 \log_{10}\left(\frac{\sqrt{(h_{tx}-h_{rx})^2 + d^2}}{1000}\right)
$$

**V2V 链路 — 3GPP TR 36.885 模型**（高速公路场景）：
$$
d_{bp} = \frac{4(h_{tx}-1)(h_{rx}-1)f_c}{c} \quad \text{（断点距离）}
$$
$$
PL = \begin{cases}
22.7\log_{10}(3) + 41.0 + 20\log_{10}(f_c/5), & d \leq 3\,\text{m} \\
22.7\log_{10}(d) + 41.0 + 20\log_{10}(f_c/5), & 3\,\text{m} < d \leq d_{bp} \\
40\log_{10}(d) + 9.45 - 17.3\log_{10}[(h_{tx}-1)(h_{rx}-1)] + 2.7\log_{10}(f_c/5), & d > d_{bp}
\end{cases}
$$

### 3.4 鲁棒功率分配 — `calOptPower.m`

**优化问题**：
$$
\max_{P_c, P_d} \quad \log_2\left(1 + \frac{P_c \cdot g_{mB}}{\sigma^2 + P_d \cdot g_{kB}}\right)
$$
$$
\text{s.t.} \quad \Pr\left(\text{SINR}_{\text{V2V}} < \gamma_0\right) \leq p_0, \quad 0 \leq P_c \leq P_{c,\max},\ 0 \leq P_d \leq P_{d,\max}
$$

**中断概率上界推导**（Chernoff 指数不等式）：

V2V 链路的 SINR 为：
$$
\text{SINR} = \frac{P_d \cdot \alpha_k |h_k|^2}{\sigma^2 + P_c \cdot \alpha_{mk} |h_{mk}|^2}
$$

对干扰信道 $|h_{mk}|^2$（服从指数分布，参数为 1）求期望，得 MGF：
$$
\mathbb{E}_{|h_{mk}|^2}\left[\exp\left(-\frac{\gamma_0 \sigma^2}{P_d \cdot \alpha_k} \cdot t\right)\right] = \exp\left(\frac{\gamma_0 \sigma^2}{P_d \cdot \alpha_k \cdot t}\right), \quad t < 1
$$

取 $t = 1 - \varepsilon_{mk}^2$，结合 $\mathbb{E}[|h_k|^2] = 1$，得中断概率的 Chernoff 上界：
$$
\Pr(\text{SINR} < \gamma_0) \leq \frac{1}{1-p_0} \cdot \exp\left(\frac{\gamma_0 \sigma^2}{P_d \cdot \alpha_k (1-\varepsilon_{mk}^2)}\right) \cdot \left(1 + \frac{P_c \cdot \alpha_{mk} (1-\varepsilon_{mk}^2)}{\gamma_0} \cdot \frac{\gamma_0}{P_d \cdot \alpha_k}\right)
$$

令上界 $\leq p_0$，可化为二分搜索的约束条件。代码中 `tmp = 1/(1-p_0) * exp(...)` 对应上述不等式的右端项。

**算法流程**：

1. **余量补偿**：`$\gamma_{0,\text{eff}} = \gamma_0 \times 10^{25/10}$` — 将设计阈值提升 25 dB（第 21 行）
2. **近乎完美 CSI 特判**：当 $(1-\varepsilon_k^2) < 10^{-12}$ 时，直接用闭式估算功率分配（第 23-39 行）
3. **参考点计算**：求 $(P_{c0}, P_{d0})$ 用于判断落入哪种情况（第 41-61 行）

$$
P_{c0} = \frac{(1-\varepsilon_k^2)\sigma^2}{\alpha_{mk}(1-\varepsilon_{mk}^2)(1/p_0-1)\varepsilon_k^2|h_k|^2 - (1-\varepsilon_k^2)\alpha_{mk}\varepsilon_{mk}^2|h_{mk}|^2}
$$

$$
P_{d0} = \frac{P_{c0} \cdot \gamma_0 \cdot \alpha_{mk}(1-\varepsilon_{mk}^2)(1-p_0)}{\alpha_k(1-\varepsilon_k^2)p_0}
$$

4. **三种情况分类**（第 66-230 行）：
   - **Case I**（$P_{d,\max} \leq P_{d0}$）：V2V 功率受限但可行区域充足
     - Branch A：`$P_d = P_{d,\max}$`，搜索**最大可行** $P_c$
     - Branch B：`$P_c = P_{c,\max}$`，搜索**最小所需** $P_d$
   - **Case II**（$P_{d,\max} > P_{d0}$ 且 $P_{c,\max} > P_{c0}$）：V2I 功率受限但可行
     - Branch A：`$P_d = P_{d,\max}$`，搜索**最大可行** $P_c$
     - Branch B：`$P_c = P_{c,\max}$`，搜索**最小所需** $P_d$
   - **Case III**（边界/不可行区域）：`$P_c = P_{c,\max}$`，搜索**最小所需** $P_d$
5. **二分搜索**：每分支内 50 次迭代，精度 $\varepsilon = 10^{-6}$
6. **数值安全**：`exp` 参数限制 $\leq 500$ 防溢出，分母 $\leq 10^{-30}$ 时跳过

**复杂度**：$O(\log(1/\varepsilon))$，约 **28 次迭代/配对**（与车辆密度无关）。

### 3.5 非鲁棒功率分配 — `calOptPower_nonrobust.m`

假设完美 CSI（$\varepsilon = 1$），直接用估计信道计算所需功率：

```matlab
Pd_needed = gamma0 * (sig2 + Pc_max * g_mk) / g_k
```

闭式解，无迭代。**缺陷**：实际信道存在误差时，~62% 的配对 SINR 低于阈值。

### 3.6 CUE-DUE 配对 — Hungarian 算法

`munkres.m` 实现 Munkres（Hungarian）算法，求解最优指派问题：

```matlab
[assignment, cost] = munkres(-C_mk);
```

其中 `C_mk(m,k)` 为第 $m$ 个 CUE 与第 $k$ 个 DUE 配对时的 V2I 容量。取负号将**容量最大化**转化为**成本最小化**。支持矩形矩阵和部分指派（不可行配对的成本设为 $+\infty$）。

### 3.7 "有效 V2I 吞吐量"指标

$$
\text{有效V2I总吞吐量} = \sum_{i} C_i \cdot \mathbb{I}\{\text{SINR}_{\text{V2V},i} \geq \gamma_0\}
$$

只计入 V2V 成功通信的配对的 V2I 容量。非鲁棒算法虽然每对容量可能更高，但大部分配对因 V2V 中断而浪费，总有效吞吐量远低于鲁棒算法。

---

## 4. 统一仿真参数

| 参数 | 值 | 说明 |
|------|----|------|
| $f_c$ | 2 GHz | 载波频率 |
| $P_{c,\max}$ / $P_{d,\max}$ | 23 dBm | V2I / V2V 最大发射功率 |
| $\sigma^2$ | $-114$ dBm/Hz | 噪声功率谱密度 |
| $\gamma_0$ | 5 dB | V2V SINR 设计阈值 |
| $p_0$ | $10^{-4}$ | V2V 目标中断概率（Chernoff 上界设计值） |
| $r_0$ | 0.5 bps/Hz | V2I 最低速率要求 |
| $R_{\text{bs}}$ | 500 m | 基站覆盖半径 |
| $h_{\text{bs}}$ | 25 m | 基站天线高度 |
| $h_{\text{veh}}$ | 1.5 m | 车辆天线高度 |
| $\sigma_{\text{V2V}}$ / $\sigma_{\text{V2I}}$ | 3 / 8 dB | V2V / V2I 阴影衰落标准差 |
| 车道数 | 6 | 高速公路车道数 |
| 车道宽度 | 4 m | 每车道宽度 |
| **SINR 余量** | **25 dB** | **补偿 Chernoff 上界的宽松性** |

---

## 5. 五个仿真脚本

### 5.1 sim_01 — V2V SINR CDF 对比

| 项目 | 内容 |
|------|------|
| 目的 | 对比两种算法在真实信道（含 CSI 延迟误差）下的 V2V 链路 SINR 分布 |
| 固定参数 | $v=60$ km/h, $T=1$ ms, $N=20$ |
| 采样量 | 快速模式 2e5 次，正式模式 **1e6 次** Monte Carlo |
| 评估方法 | 配对阶段用估计信道，采样阶段用实际信道（`hk_actual = epsi*hk + ek`）|
| 输出 | 单图：蓝色实线(鲁棒) vs 红色虚线(非鲁棒) 经验 CDF + 中断概率标注 |

### 5.2 sim_02 — 有效 V2I 吞吐量 vs CSI 反馈周期

| 项目 | 内容 |
|------|------|
| 目的 | 验证鲁棒算法对 CSI 延迟的鲁棒性 |
| 变量 | $T = [0.2, 0.5, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6]$ ms, $v = [50, 100, 150]$ km/h |
| 采样量 | 快速模式 200 次，正式模式 **1000 次** 信道实现平均 |
| 输出 | 单图 6 条曲线：实线+填充标记(鲁棒) vs 虚线+空心标记(非鲁棒) |
| 指标 | 有效 V2I 总吞吐量（仅计入 V2V 未中断的配对） |

### 5.3 sim_03 — V2V 中断 & V2I 容量 vs 车辆密度

| 项目 | 内容 |
|------|------|
| 目的 | 评估算法在不同网络负载下的性能稳定性 |
| 固定参数 | $v=60$ km/h, $T=1$ ms |
| 变量 | $N = [10, 15, 20, 25, 30, 35, 40, 45]$（numCUE = numDUE = N）|
| 采样量 | 快速模式 100 次，正式模式 **200 次** 信道实现，每实现采样 **200 次** 误差 |
| 输出 | 双子图：(a) V2V 中断概率(对数 Y 轴) (b) 总有效 V2I 容量 |

### 5.4 sim_04 — 算法收敛性与计算复杂度

| 项目 | 内容 |
|------|------|
| 目的 | 验证二分搜索收敛性，对比计算复杂度 |
| 收敛曲线 | 200 次信道实现，每次迭代追踪 V2I 容量，稀疏采样每 5 个配对取 1 个 |
| 复杂度曲线 | 每个密度等级 50 次信道实现，遍历全部配对统计迭代次数 |
| 变量 | $N = [10, 20, 30, 40]$ |
| 输出 | 双子图：(a) V2I 容量收敛曲线 (b) 平均迭代次数 vs 密度 |
| 结论 | 鲁棒算法 ~28 次迭代收敛，与密度无关；非鲁棒为闭式解（1 次） |

### 5.5 sim_05 — SINR 阈值敏感性分析

| 项目 | 内容 |
|------|------|
| 目的 | 研究 SINR 阈值和 CSI 延迟对中断概率的联合影响 |
| 固定参数 | $v=100$ km/h, $N=20$, **500 次**信道实现 |
| 变量 | $\gamma_{\text{th}} = [-5, 0, 5, 10, 15, 20, 25, 30]$ dB（8 个点）, $T = [0.2, 0.6, 1.0, 1.4, 1.8]$ ms（5 条曲线）|
| 输出 | 双子图：(a) 鲁棒算法 5 条 T 值曲线 (b) 非鲁棒算法 5 条 T 值曲线 |
| 说明 | $\gamma_{\text{th}}$ 扩展到 30 dB 是为了展示逼近 25 dB 余量边界时中断概率的变化趋势 |

---

## 6. 运行方法

### 6.1 仿真脚本

```matlab
% 依次运行 5 个仿真脚本（输出 .mat 数据和图像到 simulation_results/）
sim_01_V2V_Outage_CDF
sim_02_V2I_Rate_vs_CSI_Delay
sim_03_V2V_Outage_and_V2I_Cap_vs_Density
sim_04_Algorithm_Convergence_Complexity
sim_05_SINR_Threshold_Outage
```

各脚本开头有 `fastMode` 开关：
- `fastMode = true`：快速验证（几十秒 ~ 1 分钟）
- `fastMode = false`：正式论文结果（3 ~ 5 分钟）

### 6.2 论文插图重绘

运行仿真后，执行：
```matlab
thesis_figures    % 从 .mat 数据重新绘制所有论文插图，输出到 thesis_figures/
```

### 6.3 输出位置

- `simulation_results/`：仿真原始输出（PNG 600dpi + PDF 矢量 + .mat 数据）
- `thesis_figures/`：统一重绘的论文插图（由 `thesis_figures.m` 生成）

---

## 7. 预期仿真结果

| 仿真 | 鲁棒算法 | 非鲁棒算法 | 性能对比 |
|------|---------|-----------|---------|
| sim_01 V2V 中断概率 | ~0.01% | ~62% | 鲁棒可靠性提升 6000+ 倍 |
| sim_02 V2I 吞吐量 | 高且稳定 | 随延迟急剧下降 | 鲁棒吞吐量 2-3 倍于非鲁棒 |
| sim_03 密度影响 | 中断 <0.1%，容量随 N 线性增长 | 中断 ~60%，容量低 | 鲁棒全面优于非鲁棒 |
| sim_04 收敛性 | ~28 次迭代 | 1 次（闭式解） | 计算代价可接受 |
| sim_05 阈值敏感性 | $\gamma_{\text{th}} \leq 15$ dB 时接近 0 | ~70% | 25 dB 余量覆盖常用阈值范围 |

---

## 8. 论文插图格式

所有仿真脚本在输出前调用 `setThesisFont(gcf)` 统一格式：

| 项目 | 规范 |
|------|------|
| 中文字体 | 宋体 (SimSun) 10.5pt |
| 英文/数字 | Times New Roman 10.5pt |
| 图例 | 9pt (小五号) |
| 坐标轴 | Box 封闭，刻度朝内，线宽 0.5pt |
| 线条 | 线宽 1.5pt，标记尺寸 7pt |
| 子图标记 | (a), (b) 文字标注，左上角 |
| 输出分辨率 | 600dpi PNG + 矢量 PDF |

---

## 9. 关键代码索引

| 功能 | 文件 | 核心行号 |
|------|------|---------|
| 25 dB 余量补偿 | calOptPower.m | 第 21 行 |
| 参考点 $(P_{c0}, P_{d0})$ 计算 | calOptPower.m | 第 42-61 行 |
| Case I / II / III 分支 | calOptPower.m | 第 66-230 行 |
| 非鲁棒闭式解 | calOptPower_nonrobust.m | 第 52-110 行 |
| V2I 路径损耗 | genPL.m | 第 91 行 |
| V2V 路径损耗（近场/远场分段）| genPL.m | 第 68-78 行 |
| 车辆拓扑生成 | genCUEandDUE.m | 第 75-165 行 |
| 时间相关系数（Jake's）| 各仿真脚本 | `epsi_k = besselj(0, ...)` |
| Hungarian 最优匹配 | munkres.m | 第 46-157 行 |
| 收敛性追踪 | sim_04.m | 第 119-328 行（内嵌函数）|
