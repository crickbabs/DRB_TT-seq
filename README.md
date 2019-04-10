# DRB_TT-seq

Scripts for the analysis of TT-seq and DRB/TT-seq data.

This is a companion repository to the publication below.  

*Nascent transcriptome profiles and measurement of transcription elongation using TT-seq.*<br>
*Lea H. Gregersen<sup>1</sup> Richard Mitter<sup>2</sup> and Jesper Q. Svejstrup<sup>1</sup>.*<br>
*<sup>1</sup>Mechanisms of Transcription Laboratory, The Francis Crick Institute, 1 Midland Road, London, NW1 1AT, UK.*<br>
*<sup>2</sup>Bioinformatics and Biostatistics, The Francis Crick Institute, 1 Midland Road, London NW1 1AT, UK.*<br>

---

### align.md
This is a bash script for aligning paired Illumina sequence reads against a reference genome using STAR resulting in a sorted and indexed BAM file.

### bigwig.md
This is a bash script for generating scaled strand-specific BIGWIG files from a BAM file containing paired reads.

### metaprofiles.md
This is a bash script for generating strand-specific metagene and TSS profiles from a BAM file using "ngs.plot".

### DRB-TTseq.Rmd
This is an R markdown document describing a pipeline for calling RNA Pol II transcription wave peak positions and elongation rates from DRB/TT-seq time-series data using R.  Instructions are given for calculating wave peaks at both the single-gene and meta-gene level.
