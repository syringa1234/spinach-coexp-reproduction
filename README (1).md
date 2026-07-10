# Reproducing a Plant Saponin Gene Co-expression Network
### 复现植物三萜皂苷 SOAP 基因共表达网络 (Nature Chemical Biology, 2020)

> An end-to-end reproduction of the gene co-expression network (Fig. 1b) from
> Jozwiak et al., *"Plant terpenoid metabolism co-opts a component of the cell wall
> biosynthesis machinery"*, **Nature Chemical Biology** 16, 740–748 (2020).
>
> 本项目独立复现了该论文中菠菜 (*Spinacia oleracea*) 三萜皂苷生物合成通路的
> **基因共表达网络图 (Fig 1b)**，涵盖从原始测序数据到最终网络图的完整分析流程。

**What makes this more than a "run-the-pipeline" exercise:** the reproduction did
not match the paper on the first attempt. Rather than stopping there, I treated each
discrepancy as something to investigate. This README is organized around that
process — **what went wrong, how I found it, how I fixed it, and what I explored next.**

> 本项目的重点不在于"跑通流程"，而在于**复现结果与原文不一致时的排查、修正与延伸探索**。
> README 按"犯错 → 发现 → 修正 → 再探索"的过程组织。

---

## Pipeline / 分析流程

```
SRA reads (5 tissues, SRR11192643–647)
        │  STAR alignment (87–94% mapping)
        ▼
featureCounts / StringTie quantification
        │  merge 5 samples
        ▼
raw count matrix ──► DESeq2 normalization ──► expression matrix
        │
        ▼
Co-expression analysis
   3 baits: SOAP1 / SOAP2 / CYP716A268v2
   similarity metrics: Pearson · Spearman · Kendall (custom Python)
        │
        ▼
Network visualization in Cytoscape
   core circle + peripheral clusters + labeled bait nodes
```

| Stage | Tools |
|-------|-------|
| Alignment 比对 | STAR |
| Quantification 定量 | featureCounts, StringTie (FPKM) |
| Normalization 归一化 | DESeq2 (R) |
| Co-expression 共表达 | CoExpNetViz (Cytoscape plugin); custom Python implementation |
| Visualization 可视化 | Cytoscape |

**Baits (3 hub genes):** SOAP1 (β-amyrin synthase), SOAP2 (CYP716A268), CYP716A268v2.
These three sit at adjacent genomic loci (107620 / 107660 / 107670), consistent with
the metabolic gene cluster described in the paper.

---

## Investigation 1 — Why Were Core Genes Missing?
## 排查一：为什么关键基因掉出了网络？

### The problem / 问题

My first reproduction used **CoExpNetViz with its default percentile thresholding**.
Several pathway genes that the paper explicitly places in the network core
(SOAP3, SOAP4, SOAP5, SOAP7) **did not appear in my network at all** — they had zero edges.

初版用 **CoExpNetViz（百分位阈值法）** 复现时，多个应位于网络核心的通路基因
（SOAP3/4/5/7）**完全未出现在网络中**，与论文 Fig 1b 不符。

### The investigation / 排查过程

1. **Computed correlations directly.** SOAP5 vs the three baits gave Pearson
   r = **0.962 / 0.928 / 0.956** — all above the nominal threshold. High enough to
   belong in the network, yet absent. A contradiction.
   → 直接计算发现 SOAP5 相关性足够高，却未入网，矛盾。

2. **Checked the run log.** 41% of genes (16,235 / 39,202) were dropped for
   **near-zero variance**, and >10% of the correlation matrix was `NaN` — low-expression
   "dead genes" were polluting the computation.
   → 日志显示 41% 基因因方差近零被丢弃，相关矩阵含大量 NaN。

3. **Filtered low-expression genes and re-ran.** The `NaN` problem disappeared, but the
   target genes were **still missing** — so the dead genes were not the real cause.
   → 过滤后 NaN 消除，但目标基因仍未入网 → 另有原因。

4. **Located the root cause.** CoExpNetViz applies a **separate percentile cutoff per
   bait**, not a single absolute threshold. The target genes' correlations fell just
   below each bait's individual cutoff, so they were excluded against all three baits.
   This is a **design characteristic of the method**, amplified by the instability of
   correlation estimates at **small sample size (n = 5)**.
   → 根因：CoExpNetViz 对每个 bait 用独立百分位门槛，目标基因对每个 bait 都恰好略低于门槛而落选。
   这是方法的固有特性，叠加 5 样本下相关系数的高度不稳定。

### The fix / 修复

I switched to the method **literally described in the paper's main text** —
an **absolute Pearson r > 0.9** cutoff — to define co-expression edges directly.
After the fix, SOAP3/4/5/7 all returned correctly to the network core, matching Fig 1b.

改用论文正文字面描述的方法——**绝对阈值 Pearson r > 0.9**——直接判定边。修复后关键基因全部正确回归核心圈。

| Gene | r vs SOAP1 / SOAP2 / CYP | Result |
|------|--------------------------|--------|
| SOAP3 | 0.969 / 0.977 / 0.961 | ✅ core |
| SOAP4 | 0.974 / 0.966 / 0.982 | ✅ core |
| SOAP5 | 0.962 / 0.928 / 0.956 | ✅ core |
| SOAP7 | 0.978 / 0.949 / 0.989 | ✅ core |

> **Note / 关键认识:** The paper's Methods section actually specifies *both* steps —
> a 5/95 percentile pre-filter **followed by** a per-bait Pearson r > 0.9 filter.
> My initial attempt only reproduced the percentile step. The absolute-threshold "fix"
> was not a departure from the paper — it was **restoring the step I had missed.**
> 论文方法其实写明了两步：百分位初筛 + r>0.9 过滤。我最初只做了第一步；所谓"修正"其实是补回论文要求的第二步。

---

## Investigation 2 — Does the Expression Metric Matter?
## 探索二：表达量定义会影响结果吗？

Co-expression depends on how "expression level" is defined. I compared three definitions:
DESeq2-normalized counts, raw FPKM, and log2(FPKM + 1).

共表达依赖于"表达量"如何定义。对比三种：DESeq2 归一化 counts、FPKM 原值、log2(FPKM+1)。

### First finding — under percentile thresholding / 初步结果（百分位法下）

| Expression metric 表达量 | Threshold | Edges |
|--------------------------|-----------|-------|
| DESeq2-normalized counts | 0.915 | 1086 |
| FPKM (raw) | 0.962 | 531 |
| log2(FPKM + 1) | 0.935 | 1128 |

Under percentile thresholding, edge overlap between metrics was **alarmingly low —
as little as ~8%** (DESeq2 ∩ log2FPKM). Raw FPKM roughly halved the edge count.
初步观察：百分位法下三种定义边重合率最低仅约 8%，FPKM 原值边数腰斩。

### Re-analysis — under the corrected absolute-threshold method / 修正方法下重算

The ~8% overlap was suspicious, so I **re-ran all three metrics under the same
absolute r > 0.9 method** (Investigation 1's corrected approach), so the comparison
would be apples-to-apples.
8% 太低，可疑。于是用修正后的绝对阈值法（同口径）重算三种定义，保证可比。

| Expression metric | Filtered genes | Edges | SOAP3/4/5/7 in core? |
|-------------------|----------------|-------|----------------------|
| DESeq2 counts | 21,609 | 1096 | ✅ all 3 baits |
| FPKM (raw) | 17,811 | 900 | ✅ all 3 baits |
| log2(FPKM + 1) | 17,171 | 750 | ✅ all 3 baits |

**Pairwise edge overlap (same r > 0.9 pipeline):**

| Comparison | Overlap |
|------------|---------|
| DESeq2 ∩ FPKM | 46.6% |
| log2FPKM ∩ FPKM | 42.1% |
| DESeq2 ∩ log2FPKM | 25.3% |

### What it means / 结论

- **The core pathway genes are robust.** Under all three metrics, SOAP3/4/5/7 stay in
  the core — the key biology does not depend on the expression definition.
  关键通路基因稳健：三种定义下 SOAP3/4/5/7 都进核心圈。
- **Much of the original "8%" was the method, not the metric.** Under a consistent
  r > 0.9 pipeline, DESeq2 ∩ log2FPKM overlap rose from **8% to 25%**. The earlier
  extreme divergence was driven substantially by CoExpNetViz's percentile judgment,
  not by the expression metric itself.
  当年"8% 超低重合"很大程度是方法造成的：同口径下重合度由 8% 升至 25%。
- With only **5 samples**, correlation remains sensitive to the metric, but the effect
  is far milder once the method is held constant.
  5 样本下表达量定义仍有影响，但固定方法后差异明显收敛。

---

## Investigation 3 — Does the Similarity Metric Matter?
## 探索三：相似度度量方法会影响结果吗？

CoExpNetViz defaults to Pearson, but "similarity" can be defined many ways. I stepped
outside the tool and computed three metrics myself, then built a network from each:
**Pearson** (linear), **Spearman** (rank), **Kendall** (rank).

CoExpNetViz 默认用 Pearson，但"相似度"有多种定义。我脱离工具，自己计算三种度量并各构一个网络：Pearson（线性）、Spearman / Kendall（秩相关）。

### A fairness problem first / 先解决一个可比性问题

The three metrics are not on the same scale. With only 5 samples, **Spearman can only
take discrete values (1.0, 0.9, 0.8, 0.7 …)** — a single swapped rank drops it a whole
step. Applying Pearson's "r > 0.9" cutoff to the rank metrics would be unfair (Kendall
values run systematically lower).

三把"尺子"刻度不同：5 样本下 Spearman 只能取离散值，秩方法直接套 0.9 不公平。

**Solution:** instead of a fixed cutoff, take the **top-N edges** for each metric
(N = Pearson's edge count, ~1096). This compares *ranking structure* — "does each metric
pick the same genes as most-correlated?" — rather than incomparable absolute values.

解决方案：不用固定阈值，每种方法各取前 N 名边（N = Pearson 边数 ≈1096），比较"排名结构"而非绝对数值。

### Results / 结果

**Equivalent cutoff at top-1096:** Pearson 0.90 · Spearman 0.90 · Kendall 0.80
(the rank metrics' lower Kendall cutoff confirms the scale mismatch).

**Target-gene recovery — how many baits each key gene connects to:**

| Gene | Pearson | Spearman | Kendall |
|------|---------|----------|---------|
| SOAP3 | 3 baits | 2 baits | 2 baits |
| SOAP4 | 3 baits | 2 baits | 2 baits |
| **SOAP5** | **3 baits** | **dropped** | **dropped** |
| SOAP7 | 3 baits | 2 baits | 2 baits |

**Pairwise edge overlap:**

| Comparison | Overlap |
|------------|---------|
| Pearson ∩ Spearman | 38.7% |
| Pearson ∩ Kendall | 39.1% |
| **Spearman ∩ Kendall** | **94.3%** |

### What it means / 结论

- **The metrics split into two camps.** The two rank methods are nearly identical
  (94%), but each overlaps the linear Pearson by only ~39%. **The decisive divide is
  "linear vs rank."** Metric choice genuinely changes which gene pairs are called
  co-expressed.
  三种方法分成两派：两个秩方法几乎一致（94%），但与 Pearson 仅约 39% 重合。分野在"线性 vs 秩"。
- **SOAP5 is the most fragile key gene.** It drops out of the network under both rank
  methods. This echoes Investigation 1, where SOAP5 was also the gene most easily lost —
  a consistent signal across the whole project.
  SOAP5 最脆弱：两种秩方法下都掉出网络，与排查一中它最易丢失的现象一致。
- The rank methods also skew the per-bait balance: under Spearman/Kendall, CYP716A268v2
  collects far fewer edges (~120) than SOAP1/SOAP2 (~470–500), another artifact of coarse
  rank resolution at n = 5.
  秩方法下 bait 严重不均衡（CYP716A268v2 边数骤降），是 5 样本下秩分辨率粗的又一副作用。

---

## What These Investigations Converge On
## 三个探索共同指向什么

Three separate perturbations — changing the **thresholding method** (Inv. 1), the
**expression metric** (Inv. 2), and the **similarity metric** (Inv. 3) — all point to the
same underlying constraint:

三个独立的扰动——判定方法、表达量定义、相似度度量——都指向同一个根本约束：

> **With only 5 samples, co-expression networks are inherently unstable, and every
> methodological choice is amplified.** The paper handles this the same way: it relies on
> an absolute r > 0.9 cutoff *and* orthogonal functional validation (VIGS, heterologous
> expression), rather than trusting the small-sample correlations alone.
>
> **5 样本下共表达网络本质不稳，任何方法选择都会被放大。** 论文正是靠绝对阈值 +
> 功能验证（VIGS、异源表达）来稳住结论，而非依赖小样本相关性本身。

This is the real lesson of the reproduction: the network is a **hypothesis-generating**
tool, and its robustness must be judged against method choices and sample size — not
taken at face value.

复现的真正收获：共表达网络是**产生假设**的工具，其稳健性必须结合方法选择与样本量来判断。

---

## Skills Demonstrated / 涉及的技能与知识

- **Linux / HPC workflow** — remote CentOS server, environment management (micromamba), Xshell/Xftp.
- **RNA-seq upstream pipeline** — SRA retrieval, STAR alignment, featureCounts / StringTie quantification.
- **Normalization** — DESeq2 median-of-ratios; understanding CPM / FPKM / TPM / log transforms and when each applies.
- **Co-expression methodology** — Pearson vs Spearman vs Kendall; percentile vs absolute thresholding; fair cross-metric comparison via rank-based top-N.
- **Critical reproduction** — diagnosing a discrepancy from logs and first-principles, distinguishing method characteristics from bugs, and validating fixes against the source.
- **Scientific communication** — bilingual documentation, structured presentation of a multi-step investigation.

---

## Repository Structure / 仓库结构

```
├── scripts/
│   ├── merge_counts.sh              # merge featureCounts into a count matrix
│   ├── normalize.R                  # DESeq2 normalization
│   ├── build_fpkm_matrix.py         # build FPKM / log2FPKM matrices from StringTie
│   ├── coexpression_absolute.py     # absolute-threshold Pearson co-expression
│   ├── compare_expression_metrics.py# Investigation 2: DESeq2 / FPKM / log2FPKM
│   ├── compare_similarity_metrics.py# Investigation 3: Pearson / Spearman / Kendall
│   └── make_groups_labels.py        # node group & label tables
├── figures/
│   ├── network_reproduced.png       # corrected network (absolute r > 0.9)
│   ├── comparison_with_paper.png    # side-by-side with Fig 1b
│   └── similarity_metric_comparison.png
├── notes/
│   └── investigation.md             # full investigation log
└── README.md
```

---

## Notes / 说明

- **This is a learning reproduction** of published work, not original research.
  本项目是对已发表工作的**学习性复现**，非原创研究。
- Sequencing data are publicly available from NCBI SRA (BioProject **PRJNA609035**,
  accessions **SRR11192643–SRR11192647**).
- CoExpNetViz's percentile thresholding is a **design characteristic** of the method,
  **not a software bug**. CoExpNetViz 百分位阈值是方法设计特性，非软件缺陷。
- Gene ID / function-name mappings were derived from the paper's supplementary datasets.
- All correlation computations at n = 5 should be read as exploratory; small-sample
  instability is discussed throughout.

---

## Reference / 参考文献

Jozwiak, A., Sonawane, P.D., Panda, S. et al. *Plant terpenoid metabolism co-opts a
component of the cell wall biosynthesis machinery.* **Nat Chem Biol** 16, 740–748 (2020).
https://doi.org/10.1038/s41589-020-0541-x
