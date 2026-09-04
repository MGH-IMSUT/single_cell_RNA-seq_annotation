#!/usr/bin/env python3
"""
Ex1.txt  ->  Ex1_clone_census.tsv

Lineage-barcode clone census from a long-format (cell, BC, UMI, RelativetoMax) table.

Pipeline
  1. build a cells x barcodes matrix of RelativetoMax values
  2. cluster barcodes on 1 - Jaccard co-occurrence (average linkage, cut 0.6) -> modules
  3. merge nested subclone barcodes into their parent module -> families
  4. rule-1 exclusive ownership: each barcode belongs to the clone holding the plurality
     of the cells where it is detected; observations of a barcode in a foreign clone's
     cell are leak-in and are removed. Iterated to a fixed point.
  5. cells whose every barcode was flagged as leak-in are reassigned (not discarded) to
     the owner clone of their strongest barcode
  6. doublet detection: a cell carrying a foreign clone's barcode set at >=70% of that
     clone's expected coverage is a doublet; doublet cells are excluded from the census
  7. per-clone census: barcode count, cell count, barcodes per cell, relative-to-max,
     weakest member detection rate

Usage
  python ex1_clone_census.py Ex1.txt Ex1_clone_census.tsv
"""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

CUT = 0.6            # linkage cut on 1 - Jaccard
NEST_CONTAIN = 0.90  # child's cells that must also carry the parent barcode
NEST_RATIO = 1.5     # parent must be seen in >= this multiple of the child's cells
NEST_MINCELLS = 3    # ignore containment for barcodes seen in fewer cells than this
DOUBLET_COV = 0.70   # coverage ratio at/above which a foreign barcode set is a doublet
MAX_ITER = 50


def load_long(path):
    """Read the long table and verify RelativetoMax == UMI / max(UMI in that cell)."""
    L = pd.read_csv(path, sep="\t")
    L.columns = ["cell", "BC", "umi", "rel"]
    cellmax = L.groupby("cell").umi.transform("max")
    exact = float((np.abs(L.rel - L.umi / cellmax) < 1e-6).mean())
    if exact < 0.999:
        raise ValueError(f"RelativetoMax is not UMI/max-per-cell ({exact:.3f} of rows match)")
    return L


def build_matrix(L):
    """cells x barcodes matrix of relative-to-max values (0 = barcode absent)."""
    cells = np.array(sorted(L.cell.unique()))
    bcs = np.array(sorted(L.BC.unique()))
    ci = {c: i for i, c in enumerate(cells)}
    bi = {b: i for i, b in enumerate(bcs)}
    V = np.zeros((len(cells), len(bcs)))
    V[L.cell.map(ci).values, L.BC.map(bi).values] = L.rel.values
    return cells, bcs, V


def jaccard_modules(V, cut=CUT):
    """Cluster barcodes on 1 - Jaccard co-occurrence; return module labels and co-counts."""
    B = V > 0
    n = B.sum(0)                                   # cells per barcode
    co = B.astype(int).T @ B.astype(int)           # co-occurrence counts
    nb = V.shape[1]
    J = np.zeros((nb, nb))
    for i in range(nb):
        union = n[i] + n - co[i]
        J[i] = np.where(union > 0, co[i] / np.maximum(union, 1), 0.0)
    np.fill_diagonal(J, 1.0)
    Z = linkage(squareform(1 - J, checks=False), "average")
    return fcluster(Z, t=cut, criterion="distance"), n, co


def merge_nested(modules, n, co):
    """Merge a module into another when its barcode is contained in the other's cells.

    Barcode j is a nested subclone barcode of i when >=90% of j's cells also carry i
    and i is seen in >=1.5x as many cells overall - i.e. i was acquired first.
    """
    fam = modules.copy()
    n_merges = 0
    for i in range(len(fam)):
        for j in range(len(fam)):
            if i == j or n[j] < NEST_MINCELLS or fam[i] == fam[j]:
                continue
            if co[i, j] / n[j] >= NEST_CONTAIN and n[i] >= NEST_RATIO * n[j]:
                fam[fam == fam[j]] = fam[i]
                n_merges += 1
    return fam, n_merges


def assign_and_strip(V, fam):
    """Rule-1 exclusive ownership, iterated to a fixed point.

    Each cell is assigned to the family of its strongest barcode; each barcode is owned
    by the clone holding the plurality of its cells; observations that disagree are
    stripped. Cells that would lose every barcode are held back so the iteration can
    reach a stable state rather than emptying them.
    """
    nb = V.shape[1]
    Vcur = V.copy()
    for _ in range(MAX_ITER):
        cc = np.where((Vcur > 0).any(1), fam[np.argmax(Vcur, axis=1)], -1)
        owner = np.full(nb, -1)
        for j in range(nb):
            lab = cc[Vcur[:, j] > 0]
            lab = lab[lab >= 0]
            if len(lab):
                owner[j] = pd.Series(lab).value_counts().idxmax()
        leak = (Vcur > 0) & (owner[None, :] >= 0) & (owner[None, :] != cc[:, None]) & (cc[:, None] >= 0)
        Vnext = Vcur.copy()
        Vnext[leak] = 0
        emptied = np.where((Vcur > 0).any(1) & ~(Vnext > 0).any(1))[0]
        Vnext[emptied] = Vcur[emptied]          # hold back, resolved in the orphan step
        if np.array_equal(Vnext > 0, Vcur > 0):
            break
        Vcur = Vnext
    return Vcur, cc, owner


def resolve_orphans(Vcur, cc, owner, fam):
    """Reassign cells whose every barcode is leak-in, then strip the residual leak-in.

    Such a cell is not contaminated - its barcodes belong to a clone other than the one
    its strongest barcode's family points to. It keeps its data and moves to the owner
    clone of its strongest barcode.
    """
    leak = (Vcur > 0) & (owner[None, :] >= 0) & (owner[None, :] != cc[:, None]) & (cc[:, None] >= 0)
    clean_any = ((Vcur > 0) & ~leak).any(1)
    orphan = np.where((Vcur > 0).any(1) & ~clean_any)[0]
    cc = cc.copy()
    for c in orphan:
        cc[c] = owner[np.argmax(Vcur[c])]
    keep_row = ~np.isin(np.arange(Vcur.shape[0]), orphan)
    Vfin = Vcur.copy()
    Vfin[leak & keep_row[:, None]] = 0
    residual = (Vfin > 0) & (owner[None, :] >= 0) & (owner[None, :] != cc[:, None]) & (cc[:, None] >= 0)
    Vfin[residual] = 0
    return Vfin, cc, len(orphan)


def call_doublets(V, Vfin, cc, owner, cells):
    """A cell holding a foreign clone's barcodes at near-full coverage carries a second cell.

    Ambient contamination delivers a fragment of a clone's barcode set at low coverage;
    a doublet delivers the set as completely as a real cell of that clone would.
    """
    hit = np.where((V > 0).sum(1) > (Vfin > 0).sum(1))[0]
    events = []
    for c in hit:
        home = cc[c]
        for f in sorted(set(owner[owner >= 0])):
            if f == home:
                continue
            fb = np.where(owner == f)[0]
            seen = [j for j in fb if V[c, j] > 0]
            if not seen:
                continue
            fcells = np.where(cc == f)[0]
            expected = (Vfin[np.ix_(fcells, fb)] > 0).mean(0).mean() if len(fcells) else 0.0
            events.append({
                "cell": cells[c], "home": f"M{int(home)}", "foreign": f"M{int(f)}",
                "n_foreign_BC": len(seen),
                "cov_ratio": (len(seen) / len(fb)) / max(expected, 1e-9),
            })
    E = pd.DataFrame(events)
    doublets = sorted(set(E.loc[E.cov_ratio >= DOUBLET_COV, "cell"])) if len(E) else []
    return E, doublets


def build_census(Vfin, cc, owner, cells, bcs, doublets):
    """Per-clone table over the doublet-free cells."""
    keep = np.array([c not in set(doublets) for c in cells])
    Vc, ccc = Vfin[keep], cc[keep]
    rows = []
    for f in sorted(set(ccc[ccc >= 0])):
        cf = np.where(ccc == f)[0]
        members = [j for j in np.where(owner == f)[0] if (Vc[cf, j] > 0).any()]
        det = Vc[np.ix_(cf, members)] > 0 if members else np.zeros((len(cf), 0), bool)
        vals = Vc[np.ix_(cf, members)]
        nz = vals[vals > 0]
        rows.append({
            "clone": f"M{int(f)}",
            "n_cells": len(cf),
            "n_BC": len(members),
            "mean_BC_per_cell": round(det.sum(1).mean(), 2) if members else 0.0,
            "median_BC_per_cell": int(np.median(det.sum(1))) if members else 0,
            "median_reltomax": round(float(np.median(nz)), 3) if nz.size else np.nan,
            "mean_reltomax": round(float(nz.mean()), 3) if nz.size else np.nan,
            "min_detect_rate": round(float(det.mean(0).min()), 3) if members else np.nan,
            "BCs": ";".join(bcs[members]),
        })
    C = (pd.DataFrame(rows)
         .sort_values(["n_cells", "n_BC"], ascending=[False, False])
         .reset_index(drop=True))
    C.insert(0, "rank", np.arange(1, len(C) + 1))
    return C


def main(in_path, out_path):
    L = load_long(in_path)
    cells, bcs, V = build_matrix(L)
    modules, n, co = jaccard_modules(V)
    fam, n_nested = merge_nested(modules, n, co)
    Vcur, cc, owner = assign_and_strip(V, fam)
    Vfin, cc, n_orphan = resolve_orphans(Vcur, cc, owner, fam)
    events, doublets = call_doublets(V, Vfin, cc, owner, cells)
    census = build_census(Vfin, cc, owner, cells, bcs, doublets)
    census.to_csv(out_path, sep="\t", index=False)

    n_raw, n_kept = int((V > 0).sum()), int((Vfin > 0).sum())
    sizes = np.array([(cc == f).sum() for f in set(cc[cc >= 0])])
    p_same = float(((sizes / len(cells)) ** 2).sum())
    print(f"input              {len(L)} observations, {len(cells)} cells, {len(bcs)} barcodes")
    print(f"modules -> families {int(modules.max())} -> {len(set(fam))} ({n_nested} nested merges)")
    print(f"leak-in stripped   {n_raw - n_kept} of {n_raw} ({100*(n_raw-n_kept)/n_raw:.1f}%)")
    print(f"orphan cells       {n_orphan} reassigned, 0 discarded")
    print(f"doublets           {len(doublets)} cells "
          f"({100*len(doublets)/len(cells):.2f}%, corrected {100*len(doublets)/len(cells)/(1-p_same):.2f}%)")
    print(f"census             {len(census)} clones, {int(census.n_cells.sum())} cells, "
          f"{int(census.n_BC.sum())} barcodes, {census.n_BC.mean():.2f} BC/clone")
    print(f"wrote              {out_path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ex1.txt",
         sys.argv[2] if len(sys.argv) > 2 else "Ex1_clone_census.tsv")
