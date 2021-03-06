# DRB_TT-seq

Scripts for the analysis of TT-seq and DRB/TT-seq data.

This is a companion repository to the publication below.  

*Using TTchem-Seq to Profile Nascent Transcription and Measuring Transcript Elongation.*<br>
*Lea H. Gregersen<sup>1</sup> Richard Mitter<sup>2</sup> and Jesper Q. Svejstrup<sup>1</sup>.*<br>
*<sup>1</sup>Mechanisms of Transcription Laboratory, The Francis Crick Institute, 1 Midland Road, London, NW1 1AT, UK.*<br>
*<sup>2</sup>Bioinformatics and Biostatistics, The Francis Crick Institute, 1 Midland Road, London NW1 1AT, UK.*<br>


---

### [align.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/align.md)
This is a bash markdown document for aligning paired Illumina sequence reads against a reference genome using STAR resulting in a sorted and indexed BAM file.

### [bigwig.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig.md)
This is a bash markdown document for generating scaled strand-specific BIGWIG files from a BAM file containing paired reads.

### [bigwig_bedtools.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig_bedtools.md)
This is a bash markdown document for generating scaled strand-specific BIGWIG files from a BAM file containing paired reads.  It is meant as an alternative to **bigwig.md** to be used on large bam files when deeptools struggles.

### [metaprofiles.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/metaprofiles.md)
This is a bash markdown document for generating strand-specific metagene, TSS and TES profiles from a BAM file using "ngs.plot".

### [DRB-TTseq.md](https://github.com/crickbabs/DRB_TT-seq/blob/master/DRB_TT-seq.md)
This is a markdown document describing a pipeline for calling RNA Pol II transcription wave peak positions and elongation rates from DRB/TT-seq time-series data using R.  Instructions are given for calculating wave peaks at both the single-gene and meta-gene level.  An Rmarkdown version of the scripts is available in the [scripts](https://github.com/crickbabs/DRB_TT-seq/blob/master/scripts)
An example html output of this script is given in **[DRB-TTseq.html](https://github.com/crickbabs/DRB_TT-seq/blob/master/DRB-TTseq.html)** - view raw, save to your desktop then open with your browser in order to view it.
Users unfamiliar with R markdown are recommended to explore it using [rstudio](https://www.rstudio.com/).

### [data](https://github.com/crickbabs/DRB_TT-seq/blob/master/data/README.md)
This directory contains details of demo FASTQ data available from the NCBI's Short Read Archive (SRA).

### [scripts](https://github.com/crickbabs/DRB_TT-seq/blob/master/scripts)
This directory contains plain script versions of the markdown documents.


---

### Dependencies and requirements

Scripts were tested using the following software versions:

* SAMtools v1.3.1
* deepTools v2.5.3
* BEDTools/2.27.1
* kentUtils
* STAR 2.5.2a
* Picard v2.1.1
* R v3.5.1 running Bioconductor version 3.7
* ngsplot v2.63

Example Bash scripts are written to be executed in a linux environment.  Scripts were tested on a linux server equipped with a 8-core Intel E5-2640 Haswell CPU running at 2.6GHz and using 8 processors and 8gb RAM.  

The R script may be run on any machine able to run R v3.5.1 or higher, though for large datasets it is recommended that at least 16gb RAM be made available to the process.

---

### References

* Andrews, S. FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc (2010).

* Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21, doi:10.1093/bioinformatics/bts635 (2013).

* Hunt, S. E. et al. Ensembl variation resources. Database (Oxford) 2018, doi:10.1093/database/bay119 (2018).

* Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics. 2010 Sep 1;26(17):2204-7.

* Lawrence, M. et al. Software for computing and annotating genomic ranges. PLoS Comput Biol 9, e1003118, doi:10.1371/journal.pcbi.1003118 (2013).

* Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079, doi:10.1093/bioinformatics/btp352 (2009).

* Mammana, A. H., J. bamsignals: Extract read count signals from bam files. R package version 1.12.11 (2016).

* Quinlan AR. BEDTools: The Swiss-Army Tool for Genome Feature Analysis. CurrProtoc Bioinformatics. 2014 Sep 8;47:11.12.1-34.

* Ramirez, F., Dundar, F., Diehl, S., Gruning, B. A. & Manke, T. deepTools: a flexible platform for exploring deep-sequencing data.  Nucleic Acids Res 42, W187-191, doi:10.1093/nar/gku365 (2014).

* Shen, L., Shao, N., Liu, X. and Nestler, E. (2014) ngs.plot: Quick mining and visualization of next-generation sequencing data by integrating genomic databases, BMC Genomics, 15, 284.

* http://broadinstitute.github.io/picard/.
