#!/usr/bin/bash

## This bash script provides a means of generating scaled strand-specific BIGWIG files from a BAM file containing paired reads.

# This alternative to "bigwig.sh" uses bedtools (Quinlan, *et al* 2014) rather than deeptools to generate bedgraph files which are in turn converted to bigwig format via Jim Kent's BigWig and BigBed tools (kent, *et al*, 2010).  It is better able to deal with large bam files.

# Dependencies:<br>
#     SAMtools     http://samtools.sourceforge.net/<br>
#     bedtools     https://bedtools.readthedocs.io/en/latest/<br>
#     kentUtils    http://bioinformatics.oxfordjournals.org/content/26/17/2204.long</br>


#### Set working, temporary and results directories.
WORKDIR="/path/to/my/working_directory/"
TMPDIR="${WORKDIR}tmp/"
BIGWIGDIR="${WORKDIR}bigwig/"
mkdir -p $TMPDIR
mkdir -p $BIGWIGDIR

#### Sample information: sample name and BAM file location.
SAMPLE="WT";
BAM="${WORKDIR}${SAMPLE}.sorted.marked.bam"

#### Chromosome size information.  An unheadered, two-column tab-delimited text file: <chromosome name> <size in bases>
CHRSIZE="${WORKDIR}chrom.sizes.txt"

#### Threads - set to take advantage of multi-threading and speed things up.
THREADS=8

#### Temporaty BAM files.
BAMFOR="${TMPDIR}${SAMPLE}.fwd.bam"     # BAM file representing reads mapping to forward strand
BAMREV="${TMPDIR}${SAMPLE}.rev.bam"     # BAM file representing reads mapping to reverse strand
BAMFOR1="${TMPDIR}${SAMPLE}.fwd1.bam"
BAMFOR2="${TMPDIR}${SAMPLE}.fwd2.bam"
BAMREV1="${TMPDIR}${SAMPLE}.rev1.bam"
BAMREV2="${TMPDIR}${SAMPLE}.rev2.bam"

#### bedgraph files.
BEDGRAPH="${TMPDIR}${SAMPLE}.bedgraph"           # BEDGRAPH file representing all reads
BEDGRAPHFOR="${TMPDIR}${SAMPLE}.for.bedgraph"    # BEDGRAPH file representing reads mapping to forward strand
BEDGRAPHREV="${TMPDIR}${SAMPLE}.rev.bedgraph"    # BEDGRAPH file representing reads mapping to reverse strand

#### Sorted bedgraph files.
BEDGRAPHSORTED="${TMPDIR}${SAMPLE}.sorted.bedgraph"           # Sorted BEDGRAPH file representing all reads
BEDGRAPHFORSORTED="${TMPDIR}${SAMPLE}.for.sorted.bedgraph"    # Sorted BEDGRAPH file representing reads mapping to forward strand
BEDGRAPHREVSORTED="${TMPDIR}${SAMPLE}.rev.sorted.bedgraph"    # Sorted BEDGRAPH file representing reads mapping to reverse strand

#### bigwig files.
BIGWIG="${BIGWIGDIR}${SAMPLE}.bigwig"           # BIGWIG file representing all reads
BIGWIGFOR="${BIGWIGDIR}${SAMPLE}.for.bigwig"    # BIGWIG file representing reads mapping to forward strand
BIGWIGREV="${BIGWIGDIR}${SAMPLE}.rev.bigwig"    # BIGWIG file representing reads mapping to reverse strand

#### Scale factor.
SCALEFACTOR=1

#### Create bigwig file for all reads.
bedtools genomecov -ibam $BAM -bg -split -scale $SCALEFACTOR > $BEDGRAPH
sort -k1,1 -k2,2 $BEDGRAPH > $BEDGRAPHSORTED
bedGraphToBigWig $BEDGRAPHSORTED $CHRSIZE $BIGWIG

#### Create bigwig file for the forward strand.
samtools view -b -f 128 -F 16 --threads $THREADS $BAM > $BAMFOR1
samtools view -b -f 80  --threads $THREADS $BAM > $BAMFOR2
samtools merge --threads $THREADS -f $BAMFOR $BAMFOR1 $BAMFOR2
samtools index $BAMFOR
bedtools genomecov -ibam $BAMFOR -bg -split -strand + -scale $SCALEFACTOR > $BEDGRAPHFOR
sort -k1,1 -k2,2 $BEDGRAPHFOR > $BEDGRAPHFORSORTED
bedGraphToBigWig $BEDGRAPHFORSORTED $CHRSIZE $BIGWIGFOR

#### Create bigwig file for the reverse strand.
samtools view -b -f 144 --threads $THREADS $BAM > $BAMREV1
samtools view -b -f 64 -F 16 --threads $THREADS $BAM > $BAMREV2
samtools merge --threads $THREADS -f $BAMREV $BAMREV1 $BAMREV2
samtools index $BAMREV
bedtools genomecov -ibam $BAMREV -bg -split -strand - -scale $SCALEFACTOR > $BEDGRAPHREV
sort -k1,1 -k2,2 $BEDGRAPHREV > $BEDGRAPHREVSORTED
bedGraphToBigWig $BEDGRAPHREVSORTED $CHRSIZE $BIGWIGREV

#### Remove temporary files.
rm $BAMFOR $BAMFFOR1 $BAMFFOR2 $BAMREV $BAMREV1 $BAMREV2 $BEDGRAPH $BEDGRAPHFOR $BEDGRAPHREV $BEDGRAPHSORTED $BEDGRAPHFORSORTED $BEDGRAPHREVSORTED
