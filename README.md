# DRB_TT-seq

Scripts for the analysis of TT-seq and DRB/TT-seq data.

This is a companion repository to the publication below.  

*Using TTchem-Seq to Profile Nascent Transcription and Measuring Transcript Elongation.*<br>
*Lea H. Gregersen<sup>1</sup> Richard Mitter<sup>2</sup> and Jesper Q. Svejstrup<sup>1</sup>.*<br>
*<sup>1</sup>Mechanisms of Transcription Laboratory, The Francis Crick Institute, 1 Midland Road, London, NW1 1AT, UK.*<br>
*<sup>2</sup>Bioinformatics and Biostatistics, The Francis Crick Institute, 1 Midland Road, London NW1 1AT, UK.*<br>


Please refer to release v1.1 (Publication) for scripts matching the manscript. 

---

### [align.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/align.md)
This is a bash markdown document for aligning paired Illumina sequence reads against a reference genome using STAR resulting in a sorted and indexed BAM file.

### [bigwig.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig.md)
This is a bash markdown document for generating scaled strand-specific BIGWIG files from a BAM file containing paired reads.

### [metaprofiles.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/metaprofiles.md)
This is a bash markdown document for generating strand-specific metagene, TSS and TES profiles from a BAM file using "ngs.plot".

### [DRB-TTseq.Rmd](https://github.com/crickbabs/DRB_TT-seq/blob/master/DRB-TTseq.Rmd)
This is an R markdown document describing a pipeline for calling RNA Pol II transcription wave peak positions and elongation rates from DRB/TT-seq time-series data using R.  Instructions are given for calculating wave peaks at both the single-gene and meta-gene level.  An example html output of this script is given in **[DRB-TTseq.html](https://github.com/crickbabs/DRB_TT-seq/blob/master/DRB-TTseq.html)**

### [data](https://github.com/crickbabs/DRB_TT-seq/blob/master/data/README.md)
This directory contains details of demo FASTQ data available from the NCBI's Short Read Archive (SRA).
