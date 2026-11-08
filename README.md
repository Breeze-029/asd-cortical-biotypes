# Cortical biotypes of autism: brain-age deviation and network-diffusion epicenter analyses

Analysis code accompanying the manuscript *"Cortical biotypes of autism show distinct brain-age,
network, and molecular signatures"* (under review).

This repository contains the two analysis pipelines that produce **Figure 4** (structural brain-age
deviation) and **Figure 5** (network diffusion model, NDM, epicenter localization):

| Folder | Manuscript figure | Content |
|---|---|---|
| `01_brain_age/` | Figure 4 | Gaussian-process brain-age prediction, age-bias correction, Brain-PAD, regional correlation |
| `02_ndm/` | Figure 5 | Structural connectome construction, NDM fitting, bootstrap epicenter selection, spin-permutation test, epicenter volume deviation |

---

## 1. Contents

### `01_brain_age/`

| Script | Purpose |
|---|---|
| `tmpRunFig301GPR.m` | Single-model GPR pipeline: nested five-fold cross-validation on TDC, final model, individual-level predicted age and Brain-PAD, age-bias correction, group effect sizes. **The Brain-PAD values reported in the paper come from this script.** |
| `Fig3_02_abc_figure_generation_v2.m` | Comparison of 9 candidate model families (ridge, ridge with quadratic terms, RBF-SVR, KRR, GPR, RVR, LSBoost, bagged trees, and GPR with outlier removal), selection by external ASD correlation, bias correction, Figure 4A/4C panels, and the 68-region W-score/Brain-PAD correlation table with FDR correction. |
| `Fig3_03_b_brainPAD_region_maps.R` | Figure 4B: cortical map of the regional Brain-PAD correlation (68 Desikan-Killiany regions). |
| `Fig3_03_b_brainPAD_region_maps_HD.R` | High-resolution re-rendering of the same map for figure assembly (optional). |

### `02_ndm/`

| Script | Purpose |
|---|---|
| `Fig4_01_ndm_epicenter_fitting.m` | Builds the morphometric similarity network (MSN) from TDC, fits the NDM for every candidate seed region and diffusion time, and runs 1,000 subject-level bootstrap resamples to quantify epicenter stability. |
| `Fig4_02_epicenter_selection_ranking.m` | Ranks all 68 regions by goodness of fit and derives the candidate epicenter clusters (quantile and fixed-TopN criteria). |
| `Fig4_03a_gof_map_all_regions.R` | Figure 5A/5C: whole-cortex GOF maps. |
| `Fig4_03b_gof_map_epicenter_selected.R` | Figure 5B/5D: GOF maps in which only the selected epicenter region is coloured. |
| `Fig4_04_epicenter_volume_density_LMM.m` | Figure 5E/5F: absolute volume deviation in epicenter versus non-epicenter regions, tested with a linear mixed model (subject random intercept). |
| `Fig4_05_ndm_permutation_test.m` | 5,000 spin-permutation tests of the best-fitting seed region (preserves spatial autocorrelation). |
| `rotate_parcellation.m` | Spherical rotation utility used by `Fig4_05`. |
| `coord_lh_ct.csv`, `coord_rh_ct.csv` | Left/right hemisphere spherical coordinates of the 68 regions. Algorithmic inputs, not study data. |

---

## 2. Data (not included)

The input tables are derivatives of the **Autism Brain Imaging Data Exchange (ABIDE)**, a public,
de-identified dataset. They are not redistributed here; request them from the corresponding author
and place them in `data/`. The `data/` folder is git-ignored on purpose.

| File | Used by | Required columns |
|---|---|---|
| `ABIDE_total2.xlsx` | `01_brain_age/*.m` | `SUB_ID`, `DX_GROUP` (1 = ASD, 2 = TDC), `AGE_AT_SCAN`, and 68 volume columns `Volume_01` … `Volume_68` |
| `ABIDE_All_Metrics_Subtypes_Summaryx.xlsx` | `tmpRunFig301GPR.m`, `Fig3_02_*_v2.m`, `Fig4_01` | `SUB_ID`, `DX_GROUP`, `Volume_Subtype` (1/2 or `SUB1`/`SUB2`), `Volume_01` … `Volume_68`; `Fig4_01` additionally needs `Area_*`, `Thickness_*`, `Curve_*`, `Abscurve_*` |
| `ABIDE_All_Metrics_Subtypes_Summary_Wscore.xlsx` | `Fig4_04`, `Fig4_05` | `SUB_ID`, `DX_GROUP`, `Volume_Subtype`, `Volume_01` … `Volume_68` (W-scores) |
| `Fig4_allROI_GOF_volume.xlsx` | `Fig4_02`, `Fig4_03a`, `Fig4_03b` | three columns `ROI`, `SUB1`, `SUB2`; 68 rows (see step 2 of the NDM pipeline) |

Convention: `DX_GROUP` uses **1 = ASD, 2 = typically developing control**. The brain-age model is
trained on controls only; ASD participants are an independent external test set. Region columns
follow the Desikan-Killiany order, odd indices = left hemisphere and even indices = right hemisphere.

---

## 3. Requirements

**MATLAB** R2024a (newer releases should also work) with:
- Statistics and Machine Learning Toolbox — `fitrgp`, `fitrlinear`, `fitrsvm`, `fitrensemble`, `TreeBagger`, `corr`, `ranksum`
- Bioinformatics Toolbox — `mafdr` (Benjamini-Hochberg FDR). Without this toolbox, replace `mafdr(..., 'BHFDR', true)` by an equivalent implementation (for example `p.adjust(p, "BH")`).

**R** ≥ 4.2 with `ggseg`, `ggplot2`, `dplyr`, `readxl`, `scales`. `ggseg` downloads the
Desikan-Killiany atlas on first installation, so that step needs an internet connection.

On Apple Silicon, an Intel build of MATLAB has to be launched through Rosetta:

```bash
arch -x86_64 /Applications/MATLAB_R2024a.app/bin/matlab -batch "run('01_brain_age/tmpRunFig301GPR.m')"
```

---

## 4. Usage

```bash
git clone <this-repository>
cd asd-cortical-biotypes
# put the four input tables listed in section 2 into data/
```

All scripts resolve paths relative to their own location, so **no path editing is required** as long
as the folder layout is preserved. Results are written to `output/` (and, for the brain-age scripts,
intermediate tables are written next to the inputs in `data/`).

### Pipeline 1 — brain-age deviation (Figure 4)

The two MATLAB scripts are independent pipelines; either can be run on its own.

```bash
# 1) individual-level predicted age, Brain-PAD and group effect sizes
matlab -batch "run('01_brain_age/tmpRunFig301GPR.m')"

# 2) model comparison, Figure 4A/4C, regional correlation table (slower: 9 models x CV)
matlab -batch "run('01_brain_age/Fig3_02_abc_figure_generation_v2.m')"

# 3) Figure 4B cortical map (needs ASD_All_ROI_PAD_Correlation_Wscore.xlsx from step 2)
Rscript 01_brain_age/Fig3_03_b_brainPAD_region_maps.R
```

Run only step 1 if individual Brain-PAD values are needed, and only step 2 + 3 if the figure panels
are needed. Do not mix intermediate values between the two pipelines: they use different model-selection
procedures by design.

### Pipeline 2 — NDM epicenter (Figure 5)

```bash
# 1) MSN + NDM fitting + 1,000 bootstrap resamples
matlab -batch "run('02_ndm/Fig4_01_ndm_epicenter_fitting.m')"
#    -> output/NDM_GOF_Maps/AllROI_FullGOF_Volume_Subtype{1,2}.xlsx
#    -> output/NDM_Epicenter_Results_NoPermutation.xlsx

# 2) merge the two per-subtype GOF tables (manual step, one minute)
#    open both files, align ROI 1-68, build three columns ROI / SUB1 / SUB2,
#    save as data/Fig4_allROI_GOF_volume.xlsx

# 3) epicenter ranking and candidate clusters
matlab -batch "run('02_ndm/Fig4_02_epicenter_selection_ranking.m')"

# 4) GOF cortical maps
Rscript 02_ndm/Fig4_03a_gof_map_all_regions.R
Rscript 02_ndm/Fig4_03b_gof_map_epicenter_selected.R

# 5) epicenter versus non-epicenter volume deviation (LMM)
matlab -batch "run('02_ndm/Fig4_04_epicenter_volume_density_LMM.m')"

# 6) spin-permutation test (5,000 permutations)
matlab -batch "run('02_ndm/Fig4_05_ndm_permutation_test.m')"
```

`Fig4_05` locates its inputs by searching the script folder and `../data`, and writes
`NDM_Permutation_Results.xlsx` to the current working directory — run it from the repository root.

---

## 5. Key parameters

**Brain age**

- Training data: TDC only (`DX_GROUP == 2`); 68 Desikan-Killiany volume features; ASD used for external testing only.
- Cross-validation: nested five-fold (outer 5 / inner 5); final model refitted on all TDC. `rng(42)`.
- Final model: Gaussian process regression (`fitrgp`), squared-exponential kernel, hyperparameters by marginal likelihood.
- Age-bias correction (Cole et al., 2017): `PAD_oof = slope x age + intercept` fitted on out-of-fold TDC predictions and subtracted from all participants; observed `slope ≈ -0.577`, `intercept ≈ 9.149`.
- Regional correlations: Pearson correlation between the 68 ASD W-scores and corrected Brain-PAD; `mafdr(..., 'BHFDR', true)`.

**NDM**

- MSN: TDC only; five morphometric indices (`Volume`, `Area`, `Thickness`, `Curve`, `Abscurve`) z-scored per participant across indices, 68 x 68 Pearson correlation, group average, diagonal removed, symmetrized, negative weights set to zero, normalized.
- Graph Laplacian: `L = I - D^(-1/2) A D^(-1/2)`; diffusion `beta = 0.5`, `t = linspace(0.1, 10, 50)`.
- Epicenter threshold: `quantile 0.85` (top 15%, about 10 regions) and fixed TopN = 10; for Subtype 2 only positive GOF values are considered.
- Bootstrap: 1,000 subject-level resamples with replacement, `rng(2024)`; MSN and Laplacian held fixed. Stability = number of resamples in which the region was retained as the best seed / 1,000.
- Permutation: 5,000 spherical rotations, `rng(2024)`; the smallest attainable p value is 1/5,001 ≈ 2.0 x 10^-4.
- LMM: `AbsDev ~ RegionType + (1|Subject)`, fit with `fitlme`.

---

## 6. Expected results (sanity check)

**Brain age (after bias correction)**

| Quantity | Value |
|---|---|
| Subtype 1 Brain-PAD | n = 428, -1.93 ± 3.30 years |
| Subtype 2 Brain-PAD | n = 620, +0.88 ± 4.38 years |
| TDC Brain-PAD | n ≈ 1,163, +0.04 ± 3.57 years |
| Group comparison | one-way ANOVA F(2, 2207) = 72.20, p ≈ 4.2 x 10^-31, eta-squared = 0.061 |
| Subtype 1 vs Subtype 2 | Welch t = -11.82, p = 2.4 x 10^-30, Hedges g = -0.706 |
| Regional FDR | 61 of 68 regions survive (all negative) |
| Model performance | TDC out-of-fold r = 0.656 (MAE 4.59, RMSE 6.49); ASD external r = 0.637 (MAE 4.63, RMSE 6.77) |

**NDM**

| Quantity | Subtype 1 | Subtype 2 |
|---|---|---|
| n | 428 | 620 |
| Best seed region | ROI 18 (right medial orbitofrontal) | ROI 62 (right temporal pole) |
| Max Spearman R | 0.747 | 0.554 |
| Spin-permutation p | 0.0004 | 0.0002 |
| Bootstrap stability | 0.303 | 0.902 |
| Candidate cluster (TopN = 10) | 15, 16, 18, 19, 25, 26, 55, 56, 58, 59 | 14, 22, 24, 25, 26, 35, 43, 44, 61, 62 |
| LMM, epicenter vs non-epicenter | p = 5.9 x 10^-18, Cliff's delta = 0.074 | p = 5.0 x 10^-13, Cliff's delta = 0.041 |

Numerical agreement to about 1e-3 is expected; small differences can arise from BLAS thread counts,
MATLAB/R versions and floating-point ordering. If values differ substantially, check the subtype
label column, the `DX_GROUP` coding, the 68-column ordering (odd = left, even = right) and whether an
older version of the input table was used.

---

## 7. Notes and limitations

1. `tmpRunFig301GPR.m` and `Fig3_02_abc_figure_generation_v2.m` are two alternative brain-age pipelines with different model-selection procedures; their intermediate values are not interchangeable.
2. `Fig4_04_epicenter_volume_density_LMM.m` contains the candidate epicenter clusters as hard-coded vectors (`epi_roi_sub1`, `epi_roi_sub2`). Update them if the threshold or the input data change.
3. `Fig4_05` writes its output to the current working directory and cannot produce p values below 2.0 x 10^-4 with 5,000 permutations.
4. The R scripts resolve the repository root from the script location. Run them with `Rscript path/to/script.R`, or set the working directory to the folder containing the script; running them from an unrelated working directory will break path resolution.
5. The R scripts export through the base `pdf()` device; do not switch to `cairo_pdf()` on systems without X11/cairo support.
6. The tables used here are derivatives of public ABIDE data. Do not redistribute raw participant data, and do not use this code for clinical diagnosis.

---

## 8. Citation

If you use this code, please cite the manuscript above. The methods implemented here follow:

- Cole, J. H. et al. Predicting brain age with deep learning from raw imaging data results in a reliable and heritable biomarker. *NeuroImage* 163, 115-124 (2017).
- Le, T. T. et al. A nonlinear simulation framework supports adjusting for age when analyzing BrainAGE. *Front. Aging Neurosci.* 10, 317 (2018) — age-bias correction.
- Raj, A. et al. A network diffusion model of disease progression in dementia. *Neuron* 73, 1204-1215 (2012).
- Váša, F. et al. Adolescent tuning of association cortex in human structural brain networks. *Cereb. Cortex* 28, 281-294 (2018) — `rotate_parcellation`.
- Desikan, R. S. et al. An automated labeling system for subdividing the human cerebral cortex. *NeuroImage* 31, 968-980 (2006) — DK-68 atlas.
- Alexander-Bloch, A. et al. On testing for spatial correspondence between maps of human brain structure and function. *NeuroImage* 178, 540-551 (2018) — spin tests.
