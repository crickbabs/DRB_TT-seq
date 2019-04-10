## This scripts provides a means of generating scaled strand-specific BIGWIG files from a BAM file containing paired reads.

Adapted from:
*Ramírez, Fidel, Devon P. Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert, Andreas S. Richter, Steffen Heyne, Friederike Dündar, and Thomas  *
*deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research (2016). doi:10.1093/nar/gkw257.  *
*https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html*

Dependencies:
    SAMtools     http://samtools.sourceforge.net/
    deepTools    https://deeptools.readthedocs.io/en/develop/index.html


#### Set working, temporary and results directories
```bash
WORKDIR="/path/to/my/working_directory/"
TMPDIR="${WORKDIR}tmp/"
BIGWIGDIR="${WORKDIR}bigwig/"
mkdir -p $TMPDIR
mkdir -p $BIGWIGDIR
```


#### Sample information: sample name and BAM file location.
```bash
SAMPLE="WT";
BAM="${WORKDIR}WT.bam"
```


#### Threads - set to take advantage of multi-threading and speed things up.
```bash
THREADS=1
```


#### Temporaty BAM files
```bash
BAMFOR="${TMPDIR}${SAMPLE}.fwd.bam"     # BAM file representing reads mapping to forward strand
BAMREV="${TMPDIR}${SAMPLE}.rev.bam"     # BAM file representing reads mapping to reverse strand
BAMFFOR1="${TMPDIR}${SAMPLE}.fwd1.bam"
BAMFFOR2="${TMPDIR}${SAMPLE}.fwd2.bam"
BAMREV1="${TMPDIR}${SAMPLE}.rev1.bam"
BAMREV2="${TMPDIR}${SAMPLE}.rev2.bam"
```


#### BIGWIG files
```bash
BIGWIG="${BIGWIGDIR}${SAMPLE}.bigwig"           # BIGWIG file representing all reads
BIGWIGFOR="${BIGWIGDIR}${SAMPLE}.for.bigwig"    # BIGWIG file representing reads mapping to forward strand
BIGWIGREV="${BIGWIGDIR}${SAMPLE}.rev.bigwig"    # BIGWIG file representing reads mapping to reverse strand
```

#### Scale factor
A scale factor is used to normalise for differences in library sizes across samples.  There are many ways to generate such a factor, such as the ratio of reads between two samples of spike-ins.  A more robust scale factor may be calculated using the "estimateSizeFactors" function from the Bioconductor package DESeq2 using sample gene count information.  Setting this to 1 indicates no scaling.
```bash
SCALEFACTOR=1
```


#### Create bigwig file for all reads
```bash
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAM -o $BIGWIG
```


#### Create bigwig file for the forward strand
Get file for transcripts originating on the forward strand.
Include reads that are 2nd in a pair (128).  Exclude reads that are mapped to the reverse strand (16)
Exclude reads that are mapped to the reverse strand (16) and first in a pair (64): 64 + 16 = 80
```bash
samtools view -b -f 128 -F 16 --threads $THREADS $BAM > $BAMF1
samtools view -b -f 80  --threads $THREADS $BAM > $BAMF2
samtools merge --threads $THREADS -f $BAMF $BAMF1 $BAMF2
samtools index $BAMF
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAMF -o $BIGWIGFOR
```


#### Create bigwig file for the reverse strand
Get the file for transcripts that originated from the reverse strand:
Include reads that map to the reverse strand (128) and are second in a pair (16): 128 + 16 = 144
Include reads that are first in a pair (64), but exclude those ones that map to the reverse strand (16)
```bash
samtools view -b -f 144 --threads $THREADS $BAM > $BAMR1
samtools view -b -f 64 -F 16 --threads $THREADS $BAM > $BAMR2
samtools merge --threads $THREADS -f $BAMR $BAMR1 $BAMR2
samtools index $BAMR
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAMR -o $BIGWIGREV
```


#### Remove temporary files
```bash
rm $BAMFOR $BAMFFOR1 $BAMFFOR2 $BAMREV $BAMREV1 $BAMREV2
```
