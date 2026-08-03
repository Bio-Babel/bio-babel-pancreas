"""Shared params, palettes, marker panels, and helpers for the GSE248369 analysis."""
import numpy as np, pandas as pd, scipy.sparse as sp
import grid_py
from grid_py.renderer import CairoRenderer

SEED = 2024
DDRTREE_RANDOM_STATE = 2016   # monocle2 reduce_dimension seed

# endocrine differentiation population feeding the alpha/beta trajectory
# (delta + epsilon are distinct minor fates -> excluded so the alpha/beta bifurcation resolves cleanly)
TRAJ_POP = [
    "Neurog3-high endocrine progenitor",
    "Fev+ endocrine precursor",
    "endocrine intermediate",
    "alpha-like endocrine",
    "Pdx1+ beta-like endocrine",
    "mature beta-like endocrine",
]
ALL_ENDOCRINE = TRAJ_POP + ["delta-like endocrine", "ghrelin/epsilon-like endocrine"]

GENO_ORDER = ["control", "Brg1Δendo;Brm+/-", "DKO"]
GENO_SHORT = {"control": "control", "Brg1Δendo;Brm+/-": "dEndo", "DKO": "DKO"}
GENO_PAL = {"control": "#4daf4a", "Brg1Δendo;Brm+/-": "#377eb8", "DKO": "#e41a1c"}

ANNOT_PAL = {
    "Neurog3-high endocrine progenitor": "#6a3d9a",
    "Fev+ endocrine precursor": "#1f78b4",
    "endocrine intermediate": "#a6cee3",
    "alpha-like endocrine": "#e31a1c",
    "Pdx1+ beta-like endocrine": "#33a02c",
    "mature beta-like endocrine": "#b2df8a",
    "delta-like endocrine": "#ff7f00",
    "ghrelin/epsilon-like endocrine": "#fb9a99",
}
ARM_PAL = {"alpha": "#e31a1c", "beta": "#1b9e77", "trunk": "#888888",
           "progenitor/precursor": "#888888"}

# curated marker programs
PROG_MARKERS  = ["Neurog3", "Sox4", "Fev", "Chga", "Chgb", "Neurod1", "Pax6", "Insm1"]
ALPHA_MARKERS = ["Gcg", "Arx", "Irx1", "Irx2", "Ttr", "Pou3f4", "Peg3", "Isl1"]
BETA_MARKERS  = ["Ins1", "Ins2", "Pdx1", "Nkx6-1", "Pax4", "Iapp", "Mafa", "Slc2a2", "Ucn3", "Nnat"]
KEY_MARKERS   = ["Neurog3", "Fev", "Gcg", "Arx", "Ins2", "Pdx1", "Nkx6-1", "Mafa"]


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def gene_id(adata, sym):
    """gene_short_name -> var_name (Ensembl id); None if absent."""
    hits = adata.var_names[adata.var["gene_short_name"].astype(str).values == sym]
    return hits[0] if len(hits) else None


def present(adata, syms):
    return [s for s in syms if gene_id(adata, s) is not None]


def logcpm_col(adata, sym, layer="counts"):
    gid = gene_id(adata, sym)
    if gid is None:
        return None
    C = dense(adata.layers[layer])
    libs = C.sum(1, keepdims=True); libs[libs == 0] = 1
    j = adata.var_names.get_loc(gid)
    return np.log1p(C[:, j] / libs[:, 0] * 1e4)


def save_pheatmap(ph, filename, width=8, height=6, dpi=200):
    gt = ph.gtable if hasattr(ph, "gtable") else ph
    r = CairoRenderer(width=width, height=height, dpi=dpi, surface_type="image", bg="white")
    st = grid_py.get_state(); st.reset(); st.init_device(r)
    grid_py.grid_draw(gt); r.write_to_png(str(filename))
    return filename


def assign_arms(adata):
    """Data-driven: identify root state, alpha/beta terminal states, the alpha/beta
    branch_point, and per-cell lineage weights from the DDRTree centroid MST.
    Returns dict with keys: root_state, alpha_state, beta_state, branch_point,
    alpha_states, beta_states, and writes adata.obs['arm'] + obsm slots."""
    import igraph as ig, monocle2py as m2
    dd = adata.uns["monocle2"]["ddrtree"]
    edges = [tuple(map(int, e)) for e in np.asarray(dd["mst_edges"])]
    K = np.asarray(dd["K"]); n_cen = K.shape[1]
    cv = np.asarray(dd["closest_vertex"]).ravel().astype(int)
    G = ig.Graph(n=n_cen, edges=edges)
    deg = np.array(G.degree())
    leaves = np.where(deg == 1)[0]

    st = adata.obs["State"].astype(int).values
    pt = adata.obs["Pseudotime"].values
    gcg = logcpm_col(adata, "Gcg"); ins = logcpm_col(adata, "Ins2")
    # centroid -> dominant state, and per-centroid mean pseudotime / markers
    cdf = pd.DataFrame({"c": cv, "state": st, "pt": pt, "gcg": gcg, "ins": ins})
    cen_state = cdf.groupby("c")["state"].agg(lambda s: s.value_counts().index[0])
    cen_pt = cdf.groupby("c")["pt"].mean()
    # root centroid = leaf with lowest mean pseudotime
    leaf_pt = {int(l): cen_pt.get(l, np.inf) for l in leaves}
    root_c = min(leaf_pt, key=leaf_pt.get)
    # terminal leaves (exclude root); alpha tip = max Gcg, beta tip = max Ins
    cen_gcg = cdf.groupby("c")["gcg"].mean(); cen_ins = cdf.groupby("c")["ins"].mean()
    term = [int(l) for l in leaves if int(l) != root_c]
    alpha_c = max(term, key=lambda c: cen_gcg.get(c, -np.inf))
    beta_c  = max(term, key=lambda c: cen_ins.get(c, -np.inf))
    path_a = set(G.get_shortest_paths(root_c, to=alpha_c)[0])
    path_b = set(G.get_shortest_paths(root_c, to=beta_c)[0])
    # per-cell weights via nearest centroid; off-path cells -> closer tip
    dist_a = np.array(G.distances(source=alpha_c)[0])
    dist_b = np.array(G.distances(source=beta_c)[0])
    w_a = np.zeros(adata.n_obs); w_b = np.zeros(adata.n_obs)
    arm = np.empty(adata.n_obs, dtype=object)
    for i, c in enumerate(cv):
        ina, inb = c in path_a, c in path_b
        if ina and inb:   w_a[i] = w_b[i] = 1.0; arm[i] = "trunk"
        elif ina:         w_a[i] = 1.0;          arm[i] = "alpha"
        elif inb:         w_b[i] = 1.0;          arm[i] = "beta"
        else:  # off both main paths -> assign to closer tip
            if dist_a[c] <= dist_b[c]: w_a[i] = 1.0; arm[i] = "alpha"
            else:                      w_b[i] = 1.0; arm[i] = "beta"
    adata.obs["arm"] = pd.Categorical(arm, categories=["trunk", "alpha", "beta"])
    adata.obsm["pseudotime"] = np.column_stack([pt, pt])
    adata.obsm["cell_weights"] = np.column_stack([w_a, w_b])   # col0=alpha, col1=beta

    # which monocle2 branch_point separates alpha-state from beta-state?
    alpha_state = int(cen_state.get(alpha_c)); beta_state = int(cen_state.get(beta_c))
    root_state = int(cen_state.get(root_c))
    branch_point = None
    for bp in (1, 2, 3):
        try:
            b = m2.build_branch_cell_dataset(adata, branch_point=bp,
                                             progenitor_method="duplicate", stretch=False)
        except Exception:
            continue
        ctab = pd.crosstab(b.obs["State"].astype(int), b.obs["Branch"])
        if alpha_state in ctab.index and beta_state in ctab.index:
            # separated if alpha_state and beta_state have their mass in different branches
            abr = ctab.loc[alpha_state].idxmax(); bbr = ctab.loc[beta_state].idxmax()
            if abr != bbr:
                branch_point = bp; break
    return dict(root_state=root_state, alpha_state=alpha_state, beta_state=beta_state,
                branch_point=branch_point, alpha_c=int(alpha_c), beta_c=int(beta_c),
                n_alpha=int(w_a.sum()), n_beta=int(w_b.sum()),
                arm_counts=pd.Series(arm).value_counts().to_dict())
