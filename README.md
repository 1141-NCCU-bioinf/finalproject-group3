[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/k4IGslji)
[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-2e0aaae1b6195c2367325f4f02e2d04e9abb55f0b24a779b69b11b9e10269abc.svg)](https://classroom.github.com/online_ide?assignment_repo_id=21940683&assignment_repo_type=AssignmentRepo)

# Principles of genome folding into topologically associating domains (Szabo et al., Science 2018)

### Goal: 

Reproduce the Hi-C contact map and TAD structures of Drosophila S2R+ cells.

### Members

- 彭啟則 112703043
- 張俊喆 112703052
- 陳勁瑋 110305002

### Demo

```bash
# 1. Activate environment
conda activate hic

# 2. Plot the heatmap (Figure 1A reproduction)
hicPlotMatrix --matrix results/matrix_corrected.cool \
              --outFileName results/Szabo_Fig1A_repro.png \
              --region 2L:9935314-12973080 \
              --log1p \
              --colorMap RdYlBu_r \
              --title "Reproduced Figure 1A (S2R+ Hi-C)" \
              --perChromosome
```

### docs

- Presentation: Bioinformatics_fianl.pdf

### data (do not upload fastq file)

- Source: Drosophila melanogaster S2R+ cell line Hi-C data from NCBI SRA (Accession: `SRR5579178`).
- Format:
  - Input: Paired-end FASTQ (`.fastq.gz`)
  - Intermediate: BAM (`merged_dedup.bam`), Juicer `.hic`
  - Final: Cooler format (`.cool`) for plotting.
- Size & Processing:
  - Original dataset size: > 600 million reads (approx. 350GB uncompressed).
  - Downsampling Strategy: Due to hardware limitations (8GB RAM), we downsampled the data to 40 million reads (10% of original) using `head -n 160000000`. This volume was sufficient to visualize distinct TAD structures at 10kb resolution.

### code

- Packages used:
  - Juicer (CPU version): For core Hi-C data processing (Alignment, Chimeric handling, Deduplication).
  - BWA / Samtools: Underlying alignment and sorting tools.
  - Pigz: Used to accelerate gzip compression/decompression (replaced `gunzip` in original scripts).
  - HiCExplorer: Used for file format conversion (`.hic` to `.cool`) and high-quality visualization.
- Analysis steps:
  1. Preprocessing: Downloaded SRA data and downsampled to 40M reads using `downsample.sh`.
  1. Alignment & Building: Ran `run.sh` (customized Juicer pipeline) to map reads to dm3 genome and generate `inter_30.hic` (MAPQ > 30).
  1. Conversion: Converted Juicer's `.hic` output to `.cool` format using `hicConvertFormat`.
  1. Visualization: Plotted the TAD structure at `chr2L:9.9M-13M` using `hicPlotMatrix` with the `Fall` (RdYlBu_r) colormap.

### results

- Reproduced Part:
  - We successfully reproduced Figure 1A, showing the Topologically Associating Domains (TADs) in the region `2L:9.9M-13M`.
  - The plot clearly shows the characteristic "triangle" structures along the diagonal, consistent with the nanocompartments described in the paper.
- Improvement / Changes:
  - Optimization for Limited Resources: Successfully adapted a pipeline designed for HPC clusters to run on a standard workstation (8GB RAM) via strategic downsampling without sacrificing structural visibility.
  - Pipeline Acceleration: Replaced single-threaded `gunzip` with parallel `pigz` in the Juicer script, significantly reducing I/O bottlenecks.
  - Visualization Flexibility: Instead of relying solely on Juicebox GUI, we implemented a CLI-based plotting workflow using HiCExplorer, allowing for reproducible, scriptable figure generation.

## References

- Primary Paper: Szabo, Q., et al. "Principles of genome folding into topologically associating domains." Science 360.6385 (2018): eaar8082.
- Tools:
  - Durand, N. C., et al. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." Cell systems 3.1 (2016): 95-98.
  - Wolff, J., et al. "Galaxy HiCExplorer: a web server for reproducible Hi-C data analysis, quality control and visualization." Nucleic acids research 46.W1 (2018): W11-W16.
