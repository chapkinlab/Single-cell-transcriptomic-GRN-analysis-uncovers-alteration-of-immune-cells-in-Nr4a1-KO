#!/usr/bin/env python3
"""
plot_grn_networkx.py

Plot the two GRN adjacency matrices exported from MATLAB (build_grn_compare.m)
side by side with NetworkX, using a shared layout so node positions match
across the two conditions. Everything is a command-line argument.

  - keep the top (1-q) strongest links (|weight| >= q-th percentile, pooled)
  - shared layout (union of both filtered networks)
  - blue edges = positive weight, orange edges = negative weight
  - nodes = translucent blue discs with black gene labels on top
  - left = Nr4a1 KO (A1), right = WT (A2)

Inputs (same folder, chosen by --tag):
  A1_KO_<tag>.csv, A2_WT_<tag>.csv, gene_list_<tag>.csv
  (+ master_coords_<tag>.csv when --layout matlab)

Layouts: kamada_kawai (cleanest), spring, graphviz_neato, matlab (MATLAB coords)

Examples
--------
  python plot_grn_networkx.py --tag stem --layout kamada_kawai
  python plot_grn_networkx.py --tag stem --layout matlab
"""

import argparse
import sys
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LAYOUTS = ["kamada_kawai", "spring", "graphviz_neato", "matlab"]

POS_COLOR  = (0.0, 0.4470, 0.7410)     # blue  (positive edges)
NEG_COLOR  = (0.8500, 0.3250, 0.0980)  # orange (negative edges)

# Node color presets — deliberately NOT blue/orange, so nodes read as
# independent objects rather than part of the edge (direction) coding.
NODE_PRESETS = {
    "gray":   (0.60, 0.60, 0.60),   # neutral, recedes (default)
    "slate":  (0.35, 0.40, 0.47),   # cool dark gray
    "green":  (0.466, 0.674, 0.188),# distinct 3rd hue
    "purple": (0.494, 0.184, 0.556),
    "teal":   (0.30, 0.69, 0.68),
    "sand":   (0.85, 0.73, 0.51),   # warm neutral
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Plot WT vs Nr4a1-KO GRNs with NetworkX (shared layout).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--tag", default="stem",
                   help="dataset tag: 'stem' (Fig 4B) or 'cd8t' (Fig 6B)")
    p.add_argument("--layout", default="kamada_kawai", choices=LAYOUTS + ["all"],
                   metavar="NAME",
                   help="layout: " + ", ".join(LAYOUTS) + ", or 'all' to render every layout")
    p.add_argument("--top", type=float, default=None,
                   help="keep this fraction of strongest links (pooled); "
                        "default per tag: stem=0.25, cd8t=0.30")
    p.add_argument("--font-size", type=float, default=10, help="gene label font size")
    p.add_argument("--node-color", default="gray",
                   help="node color: a preset (" + ", ".join(NODE_PRESETS) +
                        ") or any matplotlib color / hex like '#7f7f7f'")
    p.add_argument("--node-size", type=float, default=1400, help="node disc area")
    p.add_argument("--node-alpha", type=float, default=0.30, help="node disc opacity (lower = more transparent)")
    p.add_argument("--edge-width", type=float, default=2.2, help="edge line width")
    p.add_argument("--edge-alpha", type=float, default=0.28, help="edge opacity (lower = more transparent)")
    p.add_argument("--label-top", type=int, default=0,
                   help="label only the N highest-degree nodes (0 = label all)")
    p.add_argument("--figsize", default="13x6", help="figure size in inches, WxH")
    p.add_argument("--dpi", type=int, default=300, help="output PNG resolution")
    p.add_argument("--seed", type=int, default=42, help="seed for the spring layout")
    p.add_argument("--out", default=None,
                   help="output basename (default grn_compare_networkx_<tag>_<layout>)")
    p.add_argument("--list-layouts", action="store_true",
                   help="print available layouts and exit")
    return p.parse_args(argv)


def build_graph(Af, genes):
    n = len(genes)
    G = nx.Graph()
    G.add_nodes_from(genes)
    iu = np.triu_indices(n, k=1)
    for i, j in zip(*iu):
        w = Af[i, j]
        if w != 0:
            G.add_edge(genes[i], genes[j], weight=float(w))
    return G


def compute_layout(G, kind, seed, tag):
    if kind == "kamada_kawai":
        return nx.kamada_kawai_layout(G)
    if kind == "spring":
        return nx.spring_layout(G, seed=seed,
                                k=1.5 / np.sqrt(max(G.number_of_nodes(), 1)))
    if kind == "graphviz_neato":
        return nx.nx_agraph.graphviz_layout(G, prog="neato")
    if kind == "matlab":
        c = pd.read_csv(f"master_coords_{tag}.csv", header=None).to_numpy()
        return {g: (c[i, 0], c[i, 1]) for i, g in enumerate(G.nodes())}
    raise ValueError(f"Unknown layout: {kind}")


def draw(ax, G, pos, title, args):
    ax.set_title(title, fontsize=14)
    edges = list(G.edges(data=True))
    ecolors = [POS_COLOR if d["weight"] > 0 else NEG_COLOR for _, _, d in edges]
    # edges first, wide + translucent so the network reads less noisy
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=ecolors,
                           width=args.edge_width, alpha=args.edge_alpha)
    # translucent neutral discs (no hard border) so each gene sits inside a halo
    ncolor = NODE_PRESETS.get(args.node_color, args.node_color)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[ncolor],
                           node_size=args.node_size, alpha=args.node_alpha,
                           edgecolors="none")
    # black gene labels on top
    if args.label_top and args.label_top > 0:
        deg = sorted(G.degree, key=lambda x: x[1], reverse=True)
        labels = {g: g for g, _ in deg[:args.label_top]}
    else:
        labels = {g: g for g in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax,
                            font_size=args.font_size, font_color="k")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True); s.set_color("k"); s.set_linewidth(1)


def main(argv=None):
    args = parse_args(argv)
    if args.list_layouts:
        print("Available layouts:\n  " + "\n  ".join(LAYOUTS))
        return 0

    tag = args.tag
    if args.top is None:                       # per-cell-type default
        args.top = {"stem": 0.25, "cd8t": 0.30}.get(tag, 0.25)
    A1 = pd.read_csv(f"A1_KO_{tag}.csv", header=None).to_numpy()   # Nr4a1 KO
    A2 = pd.read_csv(f"A2_WT_{tag}.csv", header=None).to_numpy()   # WT
    genes = pd.read_csv(f"gene_list_{tag}.csv", header=None).iloc[:, 0].astype(str).tolist()

    n = len(genes)
    assert A1.shape == (n, n) == A2.shape, "Matrix / gene-list size mismatch"
    np.fill_diagonal(A1, 0.0)
    np.fill_diagonal(A2, 0.0)

    pooled = np.concatenate([np.abs(A1[A1 != 0]), np.abs(A2[A2 != 0])])
    thr = np.percentile(pooled, 100 * (1 - args.top))
    A1f = np.where(np.abs(A1) >= thr, A1, 0.0)
    A2f = np.where(np.abs(A2) >= thr, A2, 0.0)

    G1 = build_graph(A1f, genes)   # KO
    G2 = build_graph(A2f, genes)   # WT
    Gm = build_graph(np.maximum(np.abs(A1f), np.abs(A2f)), genes)

    try:
        w, h = (float(v) for v in args.figsize.lower().split("x"))
    except Exception:
        w, h = 13, 6

    layouts_to_run = LAYOUTS if args.layout == "all" else [args.layout]
    print(f"tag={tag} | threshold |w| >= {thr:.4g} | "
          f"KO edges: {G1.number_of_edges()} | WT edges: {G2.number_of_edges()}")
    for lay in layouts_to_run:
        try:
            pos = compute_layout(Gm, lay, args.seed, tag)
        except Exception as e:
            print(f"[warn] layout '{lay}' failed ({e}); falling back to spring")
            pos = compute_layout(Gm, "spring", args.seed, tag)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(w, h))
        draw(ax1, G1, pos, "$\\it{Nr4a1}$ KO", args)
        draw(ax2, G2, pos, "WT", args)
        fig.patch.set_facecolor("w")
        plt.tight_layout()

        # when rendering all, always tag the file by layout so nothing overwrites
        if args.out and args.layout != "all":
            base = args.out
        else:
            base = f"grn_compare_networkx_{tag}_{lay}"
        fig.savefig(base + ".png", dpi=args.dpi, bbox_inches="tight")
        fig.savefig(base + ".pdf", bbox_inches="tight")
        fig.savefig(base + ".svg", bbox_inches="tight")
        plt.close(fig)
        print(f"  {lay:16s} -> {base}.png / .pdf / .svg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
