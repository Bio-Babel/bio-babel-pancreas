# GSE248369 — Endocrine differentiation trajectory across control / dEndo / DKO

> **Release note.** This is the preliminary report produced by the agent and is
> retained as an analysis artifact, not as the final biological interpretation.
> File names and paths below refer to the original working directory. The
> reviewed results in `followup.ipynb` and `tables/` supersede its exploratory
> genotype section. One of the 1,116 genes lacked a finite association result,
> so the statement below that all genes converged should read 1,115/1,116.
> Mechanistic suggestions in this report, including those involving *Smarca1*,
> are hypotheses rather than established conclusions.

**Analyst notebook:** `GSE248369_biobabel_analysis.ipynb`  ·  **Figures:** `outputs/figures/`  ·  **Tables:** `outputs/*.csv`
**Seed:** 2024 (Monocle2 DDRTree `random_state=2016`; tradeSeq `n_knots=6`).
**Packages (Bio-Babel):** monocle2-python 2.9.0, ddrtree-python 0.1.6, tradeSeq-python 1.13.12, ggplot2-python 4.0.2.
**Source study:** Davidson RK et al., *Diabetologia* 2024;67(10):2275–2288 (doi:10.1007/s00125-024-06211-7); data GSE248369 — E15.5 *Neurog3*-lineage (Tomato⁺) pancreatic endocrine cells, three genotypes: `control`, `Brg1Δendo;Brm+/-` (**dEndo**, endocrine BRG1 loss), `DKO` (**BRG1+BRM** double knockout). BRG1/BRM are the ATPase subunits of the SWI/SNF chromatin-remodelling complex.

---

## 1. Trajectory

From the hand-annotated object I subset the **endocrine differentiation population** (3,135 cells after dropping the `low_transcriptome_mapping_flag` cells): *Neurog3*-high progenitor (442), *Fev*⁺ precursor (819), endocrine intermediate (390), alpha-like (553), Pdx1⁺ beta-like (640), mature beta-like (291). Delta- and epsilon/ghrelin-like cells were **excluded** from the trajectory (see Limitations).

A Monocle2 **DDRTree** trajectory was built from raw counts through the canonical 9-step pipeline. The key methodological choice: ordering genes were selected **supervised** (one-vs-rest markers across the endocrine annotations) rather than by dispersion. With dispersion genes the 2-D tree is dominated by a generic maturation axis and alpha/beta do **not** separate (they share one arm); the supervised panel orients the principal graph along *fate identity* and cleanly resolves the bifurcation.

The tree was rooted automatically at the *Neurog3*-high progenitor-enriched state (State 4, lowest pseudotime, highest *Neurog3*). It forms a single main bifurcation:

```
Neurog3-high progenitor ──▶ Fev+ precursor ──▶  ┌─ ALPHA arm (Gcg⁺, State 1)
   (root, State 4)          (trunk, State 2)     └─ BETA  arm (Ins⁺, State 5)
```

Cell assignment (from the centroid MST root→tip paths): **trunk 1,213 / alpha 997 / beta 925**. A minor early side-pocket (State 3, *Fev*⁺, alpha-biased) is a secondary branch off the trunk. Figures: `ab_pseudotime.png`, `ab_annotation.png`, `ab_genotype.png`, `ab_state.png`, and per-marker overlays `ab_marker_{Neurog3,Fev,Gcg,Arx,Ins2,Pdx1}.png`.

## 2. Branch identities

Each arm carries its expected terminal program (per-arm marker-program heatmap `marker_program_heatmap.png`; z-scored, pseudotime-binned):

- **Progenitor/precursor (shared trunk):** *Neurog3*, *Sox4* high earliest; *Fev*, *Chga*, *Chgb*, *Neurod1* peak mid-trunk.
- **Alpha arm (late):** *Gcg*, *Arx*, *Irx1*, *Irx2*, *Ttr*, *Pou3f4*, *Peg3*, *Isl1*.
- **Beta arm (late):** *Ins1*, *Ins2*, *Iapp*, *Mafa*, *Slc2a2*, *Ucn3*, *Nnat*.

Quantitatively, in the last pseudotime bins the alpha program is high only in the alpha arm (mean z ≈ +0.87 vs −0.45 for the beta program) and vice-versa in the beta arm (beta program z ≈ +1.15 vs −0.47). The mapping is unambiguous: **lineage 1 = alpha, lineage 2 = beta**.

## 3. Dynamic genes

**Monocle2 BEAM** at the alpha/beta branch point flags **2,850 branch-dependent genes** (qval < 1e-4 of 16,888 tested; `beam_alpha_beta.csv`, heatmap `beam_branched_heatmap.png`). Top hits: *Wnk3, Nnat, Slc38a5, Dlg2, Ero1lb, **Mafa**, **Ttr**, **Smarca1**, Ppp1r1a, Atp2a2, Dlk1* — i.e. the canonical fate markers plus secretory/ER-stress machinery.

**tradeSeq** NB-GAMs along the two lineages (1,116 genes fit, all converged; `tradeseq_*.csv`):
- **`diffEndTest`** (different endpoint expression): *Ins1, Gcg, Ins2, Sytl4, Sec61b, C2cd4b, Nnat, Gng12, Ghrl, Dlg2, Hadh, **Smarca1*** — insulin vs glucagon dominate, as expected.
- **`patternTest`** (different trajectory shape): *Ins1, Gcg, Ins2, Nnat, Gng12, Sec61b, Gast, Atp2a2, **Smarca1**, Hadh, **Pdx1***.
- Smoothers for key markers (`tradeseq_smoothers_markers.png`) and a fitted-smoother dynamic heatmap (`tradeseq_dynamic_heatmap.png`) show *Neurog3*/*Fev* falling early, *Gcg/Arx* rising only on the alpha smoother, *Ins1/2/Mafa/Iapp* rising only on the beta smoother.

## 4. Genotype comparison

**Embryo-level allocation (the correct biological-replicate test).** Per-embryo beta fraction among committed cells, tested across genotypes:

| genotype | embryos (post-QC) | mean beta fraction |
|---|---|---|
| control | 3 | ~0.53 |
| dEndo | 4 | ~0.44 |
| DKO | 2 | ~0.46 |

Kruskal–Wallis **H = 1.71, p = 0.43** (pairwise Mann–Whitney all n.s.). → **No significant genotype difference in how cells are allocated to the alpha vs beta arms** (`allocation_per_embryo.csv`, `allocation_per_embryo.png`). River plots (`river_arm_composition.png`, `river_alpha_beta_committed.png`) show the same qualitative composition in all three genotypes (mutants trend slightly alpha-biased late, but not significantly). This **reproduces the published conclusion** that E15.5 cell allocation/trajectory is largely unperturbed.

> A naïve **cell-level** chi-square is significant (χ² = 9.98, p = 0.007) — but that treats each cell as an independent replicate (pseudoreplication) and is *not* a valid genotype test; it is reported only to illustrate the discrepancy.

**Condition-aware dynamics.** tradeSeq was refit with **genotype as the condition** (condition-specific smoothers per lineage — 6 curves/gene: lineage 1=alpha, 2=beta × control/dEndo/DKO) and `conditionTest` run per lineage. Of the **150 genes fit, 113 show genotype-dependent expression dynamics** (p < 0.05; `tradeseq_condition_test.csv`, smoothers in `tradeseq_condition_smoothers.png`). Top hits are dominated by hormone / secretory-maturation genes — ***Ppy*, Tuba1a, Cda, *Ghrl*, Mbnl2, Pde4d, Pdia3, *Pcsk1*** (and *Ins1/Ins2* among the leaders). So the *dynamics* of the secretory program along each arm differ by genotype even though *allocation* (§4) does not: not "how many cells go each way" but "does the expression trajectory along each arm shift by genotype." The smoother panels for *Ins1/Ins2/Ppy* show the control, dEndo and DKO curves separating at late pseudotime on the beta and alpha arms. **Caveat:** with DKO at 2 embryos, many hits are likely inflated by DKO sparsity; interpret primarily as **control vs dEndo** and treat the list as hypothesis-generating. That *Pcsk1* (prohormone convertase), *Ppy*, *Ghrl* and insulin genes head the list is consistent with the paper's theme that SWI/SNF loss impairs the hormone-**maturation** program more than fate allocation per se.

## 5. What the data show that the published annotations don't capture

1. **An alpha-committing "endocrine intermediate" state.** The hand annotation carries an `endocrine intermediate` population (n=390; *Arx*⁺, *Gcg*-intermediate, insulin-low) that the trajectory places firmly on the **alpha lineage** (376/390 cells map to the alpha arm) — i.e. an *Arx*⁺ pre-alpha intermediate, not a generic "intermediate." The paper's cluster scheme does not resolve this as an alpha-committed step.
2. **SWI/SNF component dynamics along the bifurcation.** ***Smarca1*** (a SWI/SNF-family ATPase) is a top BEAM branch-divergent gene and a top tradeSeq `diffEnd`/`pattern` gene — the fate decision itself is accompanied by remodeller-subunit dynamics, a direct tie-in to the study's premise that the published per-cluster DE did not surface.
3. **The trajectory the paper explicitly deferred.** The authors state they did *not* perform trajectory analysis and only inferred that trajectories "appear unchanged" from cluster proportions. This analysis supplies the actual progenitor→precursor→{alpha,beta} DDRTree and a *statistically framed* (embryo-level) allocation test, confirming no significant genotype shift — while exposing the power caveat below that pooled-cell proportions hide.
4. **A hidden replicate loss.** The entire **DKO 2** embryo (422 cells) is exactly the `low_transcriptome_mapping_flag` batch; after QC, DKO has only **2 embryos**. Pooled cell-proportion analyses mask this; any genotype claim for DKO rests on 2 replicates.

## 6. Limitations

- **DKO is underpowered:** 2 surviving embryos (~70–100 committed cells/arm). DKO smoothers/allocation estimates are fragile; genotype effects are most trustworthy for control vs dEndo.
- **Delta and epsilon fates excluded** to keep the 2-D DDRTree's alpha/beta bifurcation clean and BEAM strictly two-branch. Those lineages are not represented here (they are present in the full object).
- **2-component DDRTree** is a low-dimensional embedding; fine sub-structure and rare transitions may be compressed. Pseudotime is an ordering, not real developmental time; a single E15.5 snapshot cannot show *when* genotypes diverge.
- **Gene panels are subset for tractability** (BEAM on 16,888 detected genes; tradeSeq main fit on 1,500, condition fit on ≤400 genes selected by BEAM/variance/markers). association/pattern rankings are robust; absolute q-values depend on the panel.
- Ordering-gene selection is supervised on the existing annotations, so the trajectory's *axis* is anchored to those labels by construction (it does not invent fates, but it will not contradict them).

## 7. Follow-ups

1. **Later timepoints (E17.5, P0/P6).** The paper predicts trajectory divergence would appear later; repeat this pipeline on a post-E15.5 series to test where genotype separates.
2. **Recover DKO power** (more embryos) or model embryo as a random effect in a beta-binomial GLMM rather than the embryo-level Kruskal–Wallis used here.
3. **Validate *Smarca1* (and *Mafa/Nnat*) branch dynamics** and ask whether SWI/SNF-subunit expression itself gates the alpha/beta choice.
4. **Include delta/epsilon** via a multi-branch monocle3 / higher-dimensional trajectory to place all endocrine fates in one graph.
5. **RNA velocity / CytoTRACE** as an orthogonal check on root placement and directionality.

---

## Deliverables

| File | Contents |
|---|---|
| `GSE248369_biobabel_analysis.ipynb` | Executed, self-contained notebook (parameters + seed in cell 0) |
| `outputs/figures/ab_*.png` | DDRTree trajectory: pseudotime, State, annotation, genotype, 6 marker overlays |
| `outputs/figures/marker_program_heatmap.png` | Per-arm marker-program heatmap |
| `outputs/figures/river_*.png` | Arm-composition river plots by genotype |
| `outputs/figures/allocation_per_embryo.png` | Embryo-level allocation test |
| `outputs/figures/beam_branched_heatmap.png`, `beam_branched_markers.png` | BEAM branched-expression heatmaps |
| `outputs/figures/tradeseq_{evaluatek,smoothers_markers,dynamic_heatmap,condition_smoothers}.png` | tradeSeq diagnostics, smoothers, dynamic-gene heatmap, condition smoothers |
| `outputs/beam_alpha_beta.csv` | BEAM per-gene results |
| `outputs/tradeseq_{association,start_vs_end,diff_end,pattern,condition_test}.csv` | tradeSeq Wald-test tables |
| `outputs/allocation_per_embryo.csv`, `outputs/run_metadata.json` | per-embryo allocation; run parameters |

## Reproducibility

Global seed 2024; Monocle2 DDRTree `random_state=2016`; tradeSeq `n_knots=6`. All parameters are in the notebook's first code cell (`ANALYSIS_CORES` scales the fits). Environment: `environment/` (monocle2py 2.9.0, ddrtree-python 0.1.6, tradeSeq-python 1.13.12, ggplot2-python 4.0.2, anndata 0.12, numpy 2.5). Executed on Sherlock (Slurm, 16 cores).

## Source attribution

Study design, genotypes, and the "no cluster-proportion change at E15.5" comparison are drawn from the full text of Davidson RK et al., retrieved from **PubMed Central** (PMCID PMC11912225): *The SWI/SNF chromatin remodelling complex regulates pancreatic endocrine cell expansion and differentiation in mice in vivo*, Diabetologia 2024;67(10):2275–2288 — [10.1007/s00125-024-06211-7](https://doi.org/10.1007/s00125-024-06211-7). Data: NCBI GEO GSE248369.
