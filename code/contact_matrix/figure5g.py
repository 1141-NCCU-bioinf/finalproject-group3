import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LinearSegmentedColormap, Normalize


# ===============================
# Fixed input (Fig.5G)
# ===============================
WT_BINS      = "data/GSE99105_nm_none_40000.bins.txt"
WT_CONTACTS  = "data/GSE99105_nm_none_40000.n_contact.txt"

PH_BINS      = "data/GSE99106_nm_none_40000.bins.txt"
PH_CONTACTS  = "data/GSE99106_nm_none_40000.n_contact.txt"

START_BIN = 1723
END_BIN   = 2061

OUT_FIG = "figure/figure_5G.png"


# ===============================
# Paper-style Hi-C colormap
# ===============================
paper_colors = [
    (0/16,    "#ffffff"),
    (0.15/16, "#f5f9fd"),
    (0.3/16,  "#e6f0fa"),

    (0.6/16,  "#d6e6f5"),
    (1.0/16,  "#c6dbef"),
    (1.5/16,  "#9ecae1"),

    (2/16,    "#4292c6"),
    (3/16,    "#2171b5"),
    (4/16,    "#084594"),
    (5/16,    "#08306b"),

    (6/16,    "#5a1f0e"),
    (7/16,    "#7f2704"),
    (8/16,    "#a63603"),

    (9/16,    "#c44e03"),
    (10/16,   "#dd6b20"),
    (11/16,   "#e88d4a"),

    (12/16,   "#f0b784"),
    (13/16,   "#f6d2ad"),
    (13.5/16, "#fae6cf"),

    (14/16,   "#ffffff"),
    (14.5/16, "#d9d9d9"),
    (15/16,   "#9e9ac8"),
    (16/16,   "#253494"),
]

paper_cmap = LinearSegmentedColormap.from_list(
    "paper_hic_author",
    paper_colors,
    N=256
)


# ===============================
# Utilities
# ===============================
def load_data(bin_file, contact_file):
    bins = pd.read_csv(bin_file, sep="\t")
    contacts = pd.read_csv(contact_file, sep="\t")
    return bins, contacts


def build_hic_matrix(contacts, start_bin, end_bin):
    n = end_bin - start_bin + 1
    mat = np.full((n, n), np.nan)

    mask = (
        (contacts["cbin1"] >= start_bin) &
        (contacts["cbin1"] <= end_bin) &
        (contacts["cbin2"] >= start_bin) &
        (contacts["cbin2"] <= end_bin)
    )

    sub = contacts.loc[mask]
    print(f"Using {len(sub)} contact pairs")

    for _, row in sub.iterrows():
        i = int(row["cbin1"] - start_bin)
        j = int(row["cbin2"] - start_bin)

        obs = row["observed_count"]
        exp = row["expected_count"]

        if exp > 0:
            v = np.log2(obs / exp)
            mat[i, j] = v
            mat[j, i] = v

    return mat


# ===============================
# Main: Fig.5G (correct logic)
# ===============================
def main():
    print("Loading WT...")
    wt_bins, wt_contacts = load_data(WT_BINS, WT_CONTACTS)
    wt_mat = build_hic_matrix(wt_contacts, START_BIN, END_BIN)

    print("Loading ph505...")
    ph_bins, ph_contacts = load_data(PH_BINS, PH_CONTACTS)
    ph_mat = build_hic_matrix(ph_contacts, START_BIN, END_BIN)

    assert wt_mat.shape == ph_mat.shape
    assert wt_bins.equals(ph_bins), "WT and ph bins must be identical"

    n = wt_mat.shape[0]

    # ===============================
    # Diagonal merge (paper-style)
    # ===============================
    final_mat = np.full((n, n), np.nan)

    for i in range(n):
        for j in range(n):
            if i <= j:
                final_mat[i, j] = wt_mat[i, j]   # upper triangle = WT
            else:
                final_mat[i, j] = ph_mat[i, j]   # lower triangle = ph505

    # ===============================
    # Axis coordinates (single system)
    # ===============================
    start_bp = wt_bins.loc[wt_bins["cbin"] == START_BIN, "from.coord"].values[0]
    end_bp   = wt_bins.loc[wt_bins["cbin"] == END_BIN,   "to.coord"].values[0]

    extent = [
        start_bp / 1000, end_bp / 1000,
        start_bp / 1000, end_bp / 1000
    ]

    os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)

    norm = Normalize(vmin=0, vmax=16)

    plt.figure(figsize=(7, 7))
    im = plt.imshow(
        final_mat,
        cmap=paper_cmap,
        norm=norm,
        origin="lower",
        extent=extent,
        interpolation="spline16"
    )

    plt.colorbar(im, label="Contact enrichment (log$_2$ O/E)")

    plt.xlabel("chr3R position (kb)")
    plt.ylabel("chr3R position (kb)")
    plt.title("Fig. 5G | WT (upper triangle) vs ph$^{505}$ (lower triangle)")

    plt.tight_layout()
    plt.savefig(OUT_FIG, dpi=300)
    plt.close()

    print(f"Figure saved to {OUT_FIG}")


if __name__ == "__main__":
    main()


# python3 programs/figure5g.py