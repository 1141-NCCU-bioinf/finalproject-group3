import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from matplotlib.colors import LinearSegmentedColormap, Normalize


# ===============================
# Paper-style Hi-C colormap
# ===============================
paper_colors = [
    # --- VERY LOW: white only at ~0 ---
    (0/16,    "#ffffff"),   
    (0.15/16, "#f5f9fd"),   
    (0.3/16,  "#e6f0fa"),   

    # --- LOW: clearly blue already ---
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
# Data loading
# ===============================
def load_data(bin_file, contact_file):
    bins = pd.read_csv(bin_file, sep="\t")
    contacts = pd.read_csv(contact_file, sep="\t")
    return bins, contacts


# ===============================
# Hi-C matrix construction
# ===============================
def build_hic_matrix(contacts, start_bin, end_bin):
    """
    Build log2(O/E) Hi-C matrix for the given bin range.
    """
    n_bins = end_bin - start_bin + 1
    mat = np.full((n_bins, n_bins), np.nan)

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
            value = np.log2(obs / exp)
            mat[i, j] = value
            mat[j, i] = value

    return mat


# ===============================
# Plotting
# ===============================
def plot_and_save(mat, bins, start_bin, end_bin, output_fig):
    os.makedirs(os.path.dirname(output_fig), exist_ok=True)

    start_bp = bins.loc[bins["cbin"] == start_bin, "from.coord"].values[0]
    end_bp   = bins.loc[bins["cbin"] == end_bin,   "to.coord"].values[0]

    extent = [
        start_bp / 1000, end_bp / 1000,
        start_bp / 1000, end_bp / 1000
    ]

    norm = Normalize(vmin=0, vmax=16)

    plt.figure(figsize=(6, 6))
    im = plt.imshow(
        mat,
        cmap=paper_cmap,
        norm=norm,
        origin="lower",
        extent=extent,
        interpolation="spline16"
    )

    cbar = plt.colorbar(im)
    cbar.set_label("Contact enrichment (log$_2$ scale)")

    plt.xlabel("chr2L position (kb)")
    plt.ylabel("chr2L position (kb)")
    plt.title("Hi-C map")

    plt.tight_layout()
    plt.savefig(output_fig, dpi=300)
    plt.close()


# ===============================
# High-level pipeline
# ===============================
def run_hic_plot(bin_file, contact_file, start_bin, end_bin, output_fig):
    print("Loading data...")
    bins, contacts = load_data(bin_file, contact_file)

    print("Building Hi-C matrix...")
    mat = build_hic_matrix(contacts, start_bin, end_bin)

    print("Plotting...")
    plot_and_save(mat, bins, start_bin, end_bin, output_fig)

    print(f"Done. Figure saved to {output_fig}")


# ===============================
# CLI entry point
# ===============================
def main():
    parser = argparse.ArgumentParser(
        description="Plot Hi-C contact map (log2 O/E)"
    )

    parser.add_argument(
        "--bins",
        required=True,
        help="Bin annotation file (*.bins.txt)"
    )
    parser.add_argument(
        "--contacts",
        required=True,
        help="Hi-C contact file (*.n_contact.txt)"
    )
    parser.add_argument(
        "--start-bin",
        type=int,
        required=True,
        help="Start bin index (inclusive)"
    )
    parser.add_argument(
        "--end-bin",
        type=int,
        required=True,
        help="End bin index (inclusive)"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output figure path (e.g. figure/figure_1A.png)"
    )

    args = parser.parse_args()

    run_hic_plot(
        bin_file=args.bins,
        contact_file=args.contacts,
        start_bin=args.start_bin,
        end_bin=args.end_bin,
        output_fig=args.out
    )


if __name__ == "__main__":
    main()


# python3 programs/hic_analysis.py --bins=data/GSE99105_nm_none_40000.bins.txt --contacts=data/GSE99105_nm_none_40000.n_contact.txt --start-bin=1723 --end-bin=2061 --out=figure/figure_5A.png
# python3 programs/hic_analysis.py --bins=data/GSE99104_nm_none_20000.bins.txt --contacts=data/GSE99104_nm_none_20000.n_contact.txt --start-bin=498 --end-bin=650 --out=figure/figure_1A.png

