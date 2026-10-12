# README

This repository contains the analysis code for the following publication:

> **Navigating the rapidly evolving landscape of sequencing-based spatial transcriptomics technologies: a comprehensive multi-platform benchmark**  
> Julien Mortier, Isabelle Scheyltjens, Elina Kuznecova, Kelly Lemeire, Fleur Parmentier, Arne Defauw, Frank Vernaillen, Freya Persyn, Jolien De Block, Niels Vandamme, Stefaan Derveaux, Bert De Rybel, Diether Lambrechts, and Evelien Van Hamme. 2026.

## Code overview

| Analysis | Path to code |
|---|---|
| Preprocessing pipeline (all reads) | `<technology>/pipeline/` |
| Preprocessing pipeline (100M reads) | `<technology>/subsampling/` |
| SpatialData and overview metrics | `<technology>/harpy/` |
| Random subsampling in increments of 10M reads | `<technology>/bam/` |
| Overview and subsampling figures | `plots/bam_metrics.R`, `plots/overview_metrics.R` |
| Pairwise per-gene comparisons | `plots/correlation_plots_bam.R`, `plots/correlation_plots_bam_video.ipynb` |
| Preseq analysis | `preseq/`, `<technology>/bam/preseq_*.sh`, `plots/preseq_metrics.R` |
| Gene-set overlap and cumulative contribution curves | `plots/genes.R` |
| Gene-biotype composition | `plots/gene_biotype_summary.ipynb`, `plots/gene_biotype_plots.R` |
| Probe-multiplicity analysis | `plots/probe_count_log2fc.R` |
| Capture-area architecture | `spatial_sampling/` |
| Spatial accuracy analysis | `diffusion/` |
| General spatial plots for publication | `plots/spatial_plots.ipynb` |
| Interactive Vitessce preparation | `vitessce/prepare_vitessce.ipynb` |

## Data availability

For more information about data availability, see publication.

## Dataset identifiers

| Technology | Mouse brain | Mouse liver | Mouse lung | Human lung | Rice root |
|---|---|---|---|---|---|
| BMKMANU S1000 | exp74 | — | — | — | — |
| BMKMANU S3000 | Exp99 | Exp100 | Exp101 | Exp102 | psb |
| Curio Seeker | SC11 | SC13 | SC12 | SC14 | — |
| Nova-ST | Exp65-sampleSPC014 | Exp100-sampleSPC004 | Exp93-sampleSPC002 | Exp93-sampleSPC022 | Exp100-samplePSB |
| Stereo-Seq | spc014 | spc004 | spc002 | spc022 | psb |
| Visium (FFPE) | SCA016 | SCA014 | SCA013 | SCA015 | — |
| Visium HD (FFPE) | SCA003 | SCA002 | SCA001 | SCA010 | — |
| Visium HD (FF) | SCA019 | SCA018 | SCA017 | SCA020 | — |
| Visium HD 3′ | SCA024 | SCA021, SCA023 | — | — | FPE010 |
