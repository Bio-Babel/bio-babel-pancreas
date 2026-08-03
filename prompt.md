# Single workflow prompt

The following was the sole substantive user prompt that initiated the autonomous Claude Code analysis.

The original study is Davidson RK et al., "The SWI/SNF chromatin remodelling complex regulates pancreatic endocrine cell expansion and differentiation in mice in vivo," Diabetologia 2024;67(10):2275-2288 (doi:10.1007/s00125-024-06211-7, PMCID PMC11912225).

data/analysis_input.h5ad is a mouse embryonic pancreas single-cell dataset (GSE248369) that I QC'd and annotated by hand. Start from it, and read the object and that paper to get the design, the genotypes (control / dEndo / DKO), and the annotations.

The Python environment with all the packages is at environment/

The question I care about: how endocrine progenitors move through a precursor stage toward alpha- and beta-like cells, and whether that differentiation is perturbed across the three genotypes. Use the Bio-Babel MCP packages — monocle2-python, DDRTree-python, tradeSeq-python — and make all figures with ggplot2-python.

Analyses I'd like:

- Trajectory. Subset to the endocrine differentiation population, then build and root a Monocle2 / DDRTree trajectory from the raw counts, rooted at the progenitor end. Show the DDRTree with pseudotime, key marker genes, my cell annotation, and genotype. I'm after the progenitor → precursor → alpha/beta differentiation, so resolve the alpha and beta branches and make a per-arm marker-program heatmap showing each arm matches its expected program.

- Pseudotime composition. River plots of arm composition across pseudotime, compared between genotypes.

- Genotype allocation. At the embryo level, test whether the genotypes differ in how cells are allocated to the alpha vs beta arms.

- Gene dynamics. Run a Monocle2 BEAM analysis on the alpha/beta branch point to find genes that diverge between the two fates, with the branched-expression heatmap. Then fit tradeSeq along the alpha and beta lineages, with smoothers for the key markers and a dynamic-gene heatmap along pseudotime. Also fit a condition-aware tradeSeq with genotype as the condition and run a conditionTest to ask whether the alpha/beta dynamics shift by genotype.

Deliverables:
- One clean, runnable notebook, GSE248369_biobabel_analysis.ipynb, with parameters and the random seed recorded.
- A short report: trajectory, branch identities, dynamic genes, genotype comparison, limitations, and follow-ups — and flag anything the data show that the published annotations don't capture.
- All the figures above as workflow outputs.
