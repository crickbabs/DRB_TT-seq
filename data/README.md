Example data for running these scripts are available from the NCBI's Short Read Archive [(SRA)](https://www.ncbi.nlm.nih.gov/sra/) using the accession numbers given below.

---

# DRB_TT-seq

| accession  | description |
| ---------- | ----------- |
| SRR8112728 | DRB_10min   |
| SRR8112732 | DRB_20min   |
| SRR8112736 | DRB_30min   |
| SRR8112740 | DRB_40min   |

---

# TT-seq

| accession  | description    |
| ---------- | -------------- |
| SRR8112935 | WT1_replicate1 |
| SRR8112947 | WT1_replicate2 |
| SRR8112941 | WT2_replicate1 |
| SRR8112947 | WT2_replicate2 |

Note that details of 4 paired-end fastq files are given for the example TT-seq data.  It is recommended to align each separately.  For the purposes of visualisation in the manscript, the resulting BAM files were merged prior to further analysis.  However, this step is not necessary to achieve metaprofiles using ngs.plots and analysis of just a single sample should provide reasonable results.
