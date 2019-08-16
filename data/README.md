Example data for running these scripts are available from the NCBI's Short Read Archive [(SRA)](https://www.ncbi.nlm.nih.gov/sra/) using the accession numbers given below.

---

# DRB_TT-seq

| accession  | description |
| ---------- | ----------- |
| SRR8112728 | DRB_10min   |
| SRR8112732 | DRB_20min   |
| SRR8112736 | DRB_30min   |
| SRR8112740 | DRB_40min   |

The files above represent 4 samples generated at 4 timepoints (10, 20, 30 and 40 mniutes) after DRB release.

---

# TT-seq

| accession  | description    |
| ---------- | -------------- |
| SRR8112935 | WT1_replicate1 |
| SRR8112947 | WT1_replicate2 |
| SRR8112941 | WT2_replicate1 |
| SRR8112947 | WT2_replicate2 |

Note that details of 4 paired-end fastq files are given for the example TT-seq data.  These represent two replicates of a single wild-type biological sample, with each replicate being split across 2 sequencing runs.  It is recommended to align each separately.  For the purposes of visualisation in the manscript, the resulting BAM files were merged prior to downstream analysis.  However, this step is not necessary to achieve bigwigs/metaprofiles and analysis of just a single sample should provide reasonable results.
