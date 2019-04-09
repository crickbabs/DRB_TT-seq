#!/usr/bin/bash

############################################################################################################################################################
############################################################################################################################################################
## This scripts provides a means of generating scaled strand-specific BIGWIG files from a BAM file containing paired reads.
##
## Adapted from:
## Ramírez, Fidel, Devon P. Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert, Andreas S. Richter, Steffen Heyne, Friederike Dündar, and Thomaske
## deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research (2016). doi:10.1093/nar/gkw257.
## https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
##
## Dependencies:
## SAMtools     http://samtools.sourceforge.net/
## deepTools    https://deeptools.readthedocs.io/en/develop/index.html
############################################################################################################################################################
############################################################################################################################################################


## Software
#module load SAMtools/1.3.1-intel-2016b
#module load deepTools/2.5.3


## Directories
WORKDIR="/camp/stp/babs/working/mitterr/projects/svejstrupj/lea.gregersen/SCAF.methods_paper/work/"
TMPDIR="${WORKDIR}tmp/"
BIGWIGDIR="${WORKDIR}bigwig/"
mkdir -p $TMPDIR
mkdir -p $BIGWIGDIR


## Sample information
SAMPLE="WT";
BAM="${WORKDIR}WT.bam"


## Threads
THREADS=8


## Temporaty BAM files
BAMFOR="${TMPDIR}${SAMPLE}.fwd.bam"     # BAM file representing reads mapping to forward strand
BAMREV="${TMPDIR}${SAMPLE}.rev.bam"     # BAM file representing reads mapping to reverse strand
BAMFFOR1="${TMPDIR}${SAMPLE}.fwd1.bam"
BAMFFOR2="${TMPDIR}${SAMPLE}.fwd2.bam"
BAMREV1="${TMPDIR}${SAMPLE}.rev1.bam"
BAMREV2="${TMPDIR}${SAMPLE}.rev2.bam"


## BIGWIG files
BIGWIG="${BIGWIGDIR}${SAMPLE}.bigwig"           # BIGWIG file representing all reads
BIGWIGFOR="${BIGWIGDIR}${SAMPLE}.for.bigwig"    # BIGWIG file representing reads mapping to forward strand
BIGWIGREV="${BIGWIGDIR}${SAMPLE}.rev.bigwig"    # BIGWIG file representing reads mapping to reverse strand


## Scale factor.  Used to scale the data in order to normalise across samples.  For example, these may be derived from DESeq2's estimateSizeFactors function.  Setting this to 1 indicates no scaling.
SCALEFACTOR=1


## Create bigwig file for all reads
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAM -o $BIGWIG


## Get file for transcripts originating on the forward strand.
## Include reads that are 2nd in a pair (128).  Exclude reads that are mapped to the reverse strand (16)
## Exclude reads that are mapped to the reverse strand (16) and first in a pair (64): 64 + 16 = 80
samtools view -b -f 128 -F 16 --threads $THREADS $BAM > $BAMF1
samtools view -b -f 80  --threads $THREADS $BAM > $BAMF2
samtools merge --threads $THREADS -f $BAMF $BAMF1 $BAMF2
samtools index $BAMF
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAMF -o $BIGWIGFOR


## Get the file for transcripts that originated from the reverse strand:
## Include reads that map to the reverse strand (128) and are second in a pair (16): 128 + 16 = 144
## Include reads that are first in a pair (64), but exclude those ones that map to the reverse strand (16)
samtools view -b -f 144 --threads $THREADS $BAM > $BAMR1
samtools view -b -f 64 -F 16 --threads $THREADS $BAM > $BAMR2
samtools merge --threads $THREADS -f $BAMR $BAMR1 $BAMR2
samtools index $BAMR
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAMR -o $BIGWIGREV


## Remove temporary files
rm $BAMFOR $BAMFFOR1 $BAMFFOR2 $BAMREV $BAMREV1 $BAMREV2
