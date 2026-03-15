# Multi-Response Phylogenetic Mixed Models (MR-PMM)

This repository accompanies the article **"Multi-Response Phylogenetic Mixed Models: Concepts and Application"** by Ben Halliwell, Barbara R. Holland, and Luke A. Yates.

Paper link: https://onlinelibrary.wiley.com/doi/full/10.1111/brv.70001

- Tutorial 1 (model specification, diagnostics, validation, inference, and simulation):  
  https://Benjamin-Halliwell.github.io/MR-PMM/MR-PMM_tutorial.html
- Tutorial 2 (Eucalyptus trait case study using AusTraits):  
  https://Benjamin-Halliwell.github.io/MR-PMM/MR-PMM_euc_example_analysis.html
- Additional model validation demo:  
  https://Benjamin-Halliwell.github.io/MR-PMM/modelValidation/PMM_validation_md.html

---

## 中文导读：MR-PMM 是什么？

MR-PMM（多响应系统发育混合模型）是一类专门面向生态与进化数据的统计模型。它同时建模多个性状（例如叶面积、木材密度、种子质量），并把“物种间亲缘关系（系统发育树）”纳入模型结构。

简单说，它回答三类问题：

1. **这个性状有多少“谱系惯性”**（phylogenetic signal）？
2. **两个性状相关，是因为共同祖先，还是因为当代生态过程？**
3. **多个性状之间是否存在稳健的协同/权衡关系？**

---

## 如何用 MR-PMM 做生态学分析（实践路线）

### 1) 准备三类数据

- **性状矩阵**：行为物种、列为多个性状（可连续、二分类等，视模型设定）。
- **协变量**：气候、土壤、生活型等固定效应。
- **系统发育树**：物种名称需与性状数据一致。

### 2) 拆分“相关性来源”

MR-PMM 会把性状协方差拆成至少两层：

- **系统发育层（phylogenetic level）**：由共享祖先导致的“深层相关”。
- **残差/个体层（independent/residual level）**：与亲缘关系无关、更多反映现时生态过程。

这一步对于生态解释很关键：

- 若相关主要出现在系统发育层，意味着性状关联可能长期保守。
- 若相关主要出现在残差层，说明当前环境过滤或功能策略更重要。

### 3) 读结果时看什么

- **每个性状的 phylogenetic signal**（通常用 \(h^2\) 或等价比例解释）。
- **系统发育相关 vs 非系统发育相关**（同一对性状分别报告）。
- **不确定性区间**（后验分布、可信区间），避免只看点估计。

### 4) 生态学解释模板（可直接套用）

- “性状 A 与性状 B 在系统发育层相关显著，而残差层弱，表明该权衡主要源于进化历史保守性。”
- “性状 C 的系统发育信号低，但残差层与环境变量强相关，说明其可塑性高、对局地环境响应更快。”

---

## 代码结构（给新读者的地图）

- `MR-PMM_tutorial.Rmd`：完整教学流程（从模拟到拟合与解释）。
- `MR-PMM_tutorial_functions.R`：教程中重复使用的核心函数（模拟、相关分解）。
- `MR-PMM_euc_example_analysis.Rmd`：桉树（Eucalyptus）真实案例分析。
- `MR-PMM_euc_example_functions.R`：案例分析的工具函数。
- `modelValidation/`：模型验证与诊断示例。

---

## 给大众的科普式理解（严谨但不枯燥）

把生态系统想成“家族企业联盟”：

- 物种像不同家族分支，亲缘树就是族谱；
- 性状像企业的经营策略（快周转 vs 重资产）；
- 你观测到两个策略一起出现，不一定是“今天学出来的”，也可能是“祖上传下来的”。

MR-PMM 的价值在于：

- **不把“祖传”误当“环境效应”**，
- 也不把“环境塑造”误当“进化宿命”。

因此，它特别适合用于：

- 功能性状协同/权衡研究，
- 气候适应与谱系保守性并存问题，
- 跨物种比较中“相关≠因果”的谨慎拆解。

---

## 建议阅读顺序

1. 先看 `MR-PMM_tutorial.html`，理解模型思想与术语。  
2. 再看 `MR-PMM_euc_example_analysis.html`，理解真实数据流程。  
3. 最后结合 `MR-PMM_tutorial_functions.R` 的注释，掌握可复用函数。
