# 卫星定位干扰源检测与定位研究

> 基于GNSS-R（GNSS反射测量）技术，利用CYGNSS低轨卫星星座数据，
> 研究GPS L1频段地面射频干扰源（RFI）的大范围检测与定位方法。
> 兼顾压制干扰（Jamming）与欺骗干扰（Spoofing）两类场景。

---

## 当前进展

- [x] 课题背景调研（干扰建模、低轨无源定位综述）
- [x] 核心文献阅读（CYGNSS噪声特性、RFI检测方法）
- [ ] 复现噪声底层异常检测（Level 1标准产品全球热点图）
- [ ] 复现峰度检验算法（Kurtosis, Δκ > 0.1阈值）
- [ ] 复现跨频率比较算法（Cross-frequency, Δζ）
- [ ] 探索大范围干扰发现方法
- [ ] Spoofing欺骗干扰检测方法研究

---

## 技术背景

### 核心原理：GNSS-R
GNSS信号被地面反射后由低轨卫星接收，形成**延迟多普勒图（DDM）**。
正常情况下DDM负延迟区域为纯噪声底层；地面RFI源会导致噪声底层异常抬升（最高可达8dB）。

### 数据平台：CYGNSS
- 8颗低轨小卫星星座，覆盖 ±38° 纬度
- 主要使用 **Level 1 标准产品**（DDM噪声底层统计量）
- 特殊模式：**Raw I/F原始中频数据**（约500条轨迹，每条≤60s）

### 检测方法体系

| 方法 | 原理 | 阈值 | 特点 |
|------|------|------|------|
| 噪声底层异常 | 全球噪声底层均值图，热点即干扰 | 相对均值 >几dB | 宏观普查，无需Raw数据 |
| 峰度检验（Kurtosis） | 正常信号κ=3，RFI导致偏离 | Δκ > 0.1 | 检出率~10%，适合强干扰 |
| 跨频率比较（Cross-freq） | 与L1 C/A标准谱形对比 | Δζ = 0.05~0.075 | 能检出低幅度干扰，与Kurtosis互补 |

> 两种Raw I/F算法联合使用：75.70%的时段两者一致"无RFI"，3.20%的时段同时检出RFI。

### 干扰类型
- **Jamming（压制干扰）**：地面发射机功率覆盖GPS L1频段，使接收机无法正常工作
- **Spoofing（欺骗干扰）**：伪造GNSS信号，使接收机输出错误位置/时间信息

---

## 仓库结构

```
satellite-interference-study/
├── README.md                        # 本文件（项目总览）
├── simulation/                      # 仿真模块
│   ├── code/                        # MATLAB .m 文件
│   ├── results/                     # 仿真结果（图片、数据）
│   └── notes/                       # 仿真思路记录
├── research/                        # 课题调研模块
│   ├── docs/                        # Word报告 .docx
│   ├── slides/                      # 演示文稿 .pptx
│   ├── scripts/                     # Python脚本 .py
│   └── literature/                  # 调研相关文献整理
├── literature/                      # 文献库笔记（PDF存本地/iCloud）
│   └── notes.md
└── .gitignore
```

---

## 关键文献索引

> PDF原文存于本地/iCloud「王峰老师-卫星定位干扰文献库」，此处记核心结论。

| 文献 | 核心结论 |
|------|---------|
| Al-Khaldi et al. (2025) *Detection and Analysis of GPS L1-Band RFI Using Spaceborne GNSS-R* | 用CYGNSS Raw I/F数据，联合Kurtosis+Cross-frequency算法检测RFI，约25%轨迹含明显干扰 |
| Gleason et al. (2020) *Characterizing Background Signals and Noise in Spaceborne GNSS Reflection Ocean Observations* | CYGNSS噪声底层受SBAS卫星反射影响显著（WAAS），但对Level 1定标影响可忽略 |
| 面向GNSS干扰监测的低轨星载无源定位技术综述与展望 | 低轨星座对地面干扰源无源定位的技术路线综述，涵盖TDOA/FDOA等方法 |
| 干扰建模.docx | 干扰信号建模方法，包含Jamming与Spoofing场景 |

---

## 关键术语

| 术语 | 说明 |
|------|------|
| GNSS-R | GNSS反射测量，利用反射信号进行遥感 |
| CYGNSS | 飓风全球导航卫星系统，8颗低轨小卫星，NASA |
| DDM | 延迟多普勒图，GNSS-R的基本观测量 |
| RFI | 射频干扰，Radio Frequency Interference |
| Jamming | 压制干扰，通过强功率信号淹没正常GNSS信号 |
| Spoofing | 欺骗干扰，伪造GNSS信号使接收机输出错误结果 |
| Kurtosis (κ) | 峰度，用于检测信号统计分布偏离高斯分布的程度 |
| Raw I/F | 原始中频数据，相关前的I/Q采样数据 |
| TDOA/FDOA | 到达时间差/到达频率差，无源定位常用方法 |
| SBAS/WAAS | 星基增强系统/广域增强系统，会影响CYGNSS噪声底层 |
| Level 1 | CYGNSS一级数据产品，含DDM及噪声底层统计量 |

---

## 待解决问题

- [ ] CYGNSS Level 1数据获取与预处理流程
- [ ] 噪声底层全球地图复现（对应Al-Khaldi 2025 Fig.1）
- [ ] Raw I/F数据获取途径（NASA CSDA项目）

---

## 导师
王峰老师 | 课题方向：卫星定位干扰源检测与定位
