# Papers discussing diffusion in spatial transcriptomics

## Literature

- [1] [Systematic comparison of sequencing-based spatial transcriptomic methods](https://www.nature.com/articles/s41592-024-02325-3?error=cookies_not_supported&code=5cc4c3d5-ac1a-405c-b70a-8ff096e9d1a8) See section "molecule lateral diffusion".
- [2] [Nova-ST: Nano-patterned ultra-dense platform for spatial transcriptomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC11384075/) See figure 3G.
- [3] [Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays](https://www.sciencedirect.com/science/article/pii/S0092867422003993) Stereoseq paper. Three approaches to quantize diffusion: (1) non-coding RNA (with Malat1 and Neat1), (2) mitochondriaal RNA, (3) genes specific for some cell types (Vip in Vip+ interneurons) can be used measure distance of transcripts to centroid of cells.
- [4] [Slide-seq: A Scalable Technology for Measuring Genome-Wide Expression at High Spatial Resolution](https://pmc.ncbi.nlm.nih.gov/articles/PMC6927209/)
- [5] [Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8606189/#F2) See figure 2B. Plots of the spatial expression of Hpca and Camk2a across CA1 soma and dendrites in the hippocampus.

About diffusion but not directly applicable:

- [6] [Why is no one talking about lateral diffusion in spatial sequencing?](https://www.linkedin.com/pulse/why-one-talking-lateral-diffusion-spatial-sequencing-brant-bassam-evnjf/) general comment
- [7] [From Transcripts to Cells: Dissecting Sensitivity, Signal Contamination, and Specificity in Xenium Spatial Transcriptomics](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1.full) about Xenium, inaccurate segmentation vs diffusion, SPLIT R package, using RCTD deconvolution
- [8] [Spatial genomics enables multi-modal study of clonal heterogeneity in tissues](https://pmc.ncbi.nlm.nih.gov/articles/PMC9301586/) Applicable ideas?

Miscellaneous:

- [9] [Oit3, a promising hallmark gene for targeting liver sinusoidal endothelial cells](https://www.nature.com/articles/s41392-023-01621-2) Oit3 is marker gene for LSECs? Do the LSECs form a single layer of cells? Can this be used to measure diffusion of Oit3 from the center of the nucleus of these cells, or the center of this cell line? Is it specific to LSECs, or also expressed in hepatocytes?
- [10] [xxx seq-scope xxx](https://pmc.ncbi.nlm.nih.gov/articles/PMC8238917/) See Fig 3H, genes enriched in nucleus (e.g. Malat1, Neat1 but also Mlxipl; study in liver)

## Useful genes

This leads to the following potentially useful genes for characterizing diffusion:

- Malat1 aka Neat2, and Neat1. Long non-coding RNA, expressed only in the nucleus. See [2, 3].
- Ccn2 (brain, cortical layer 6b, see [2])
- Hpca and Camk2a (brain, see [5])
- Ptgds (brain, see [1])
- Vip (brain, in Vip+ interneurons, see [3])
- Oit3 (??, liver, see [9])
