## This bash script provides a means of generating strand-specific metagene, TSS and TES profiles from a BAM file using "ngs.plot".

*Of course, this will only work if your libraries were created in a strand-specific fashion.*

Shen, L.*, Shao, N., Liu, X. and Nestler, E. (2014).<br>
ngs.plot: Quick mining and visualization of next-generation sequencing data by integrating genomic databases, BMC Genomics, 15, 284.<br>
https://github.com/shenlab-sinai/ngsplot<br>

Dependencies:<br>
SAMtools     http://samtools.sourceforge.net/<br>
ngs.plot     https://github.com/shenlab-sinai/ngsplot<br>

---


#### Set working, temporary and results directories.
```bash
WORKDIR="/path/to/my/working_directory/"
TMPDIR="${WORKDIR}tmp/"
PROFDIR="${WORKDIR}metaprofiles/"
mkdir -p $TMPDIR
mkdir -p $PROFDIR
```


#### Sample information: sample name and BAM file location.
```bash
SAMPLE="WT";
BAM="${WORKDIR}WT.bam"
MATE1="${WORKDIR}WT.mate1.bam"
```


#### Threads - set to take advantage of multi-threading and speed things up.
```bash
THREADS=1
```

#### Restict to first mate reads.
If the BAM file contains paired reads, create a new version containing only the first mate reads.<br>
This step may be skipped if your reads are not paired.<br>
```bash
samtools view --threads $THREADS -h -b -f 64 $BAM -o $MATE1
samtools index $MATE1
```


#### BAM header compliance.
It might be necessary to re-header the BAM file so that the chromosome names match those in the ngs.plot database, e.g. standard chromosomes preceeded with "chr" for hg38.  This is most easily achieved using SAMtools "reheader" function.<br>
For human hg38 Ensembl alignments the following steps should do the job.<br>
This step may be skipped if your header is already compliant.<br>
```bash
MATE1REHEADER="${WORKDIR}WT.mate1.reheader.bam"
samtools view --threads $THREADS -H ${MATE1} | sed -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - ${MATE1} > ${MATE1REHEADER}
samtools index ${MATE1REHEADER}
```


#### Run ngs.plot.
Finally, use ngs.plot to create sense and anti-sense profiles for your regions of interest using the correct BAM file.
Note that if your libraries were produced such that mate2 represents the forward strand the sense/anti-sense profiles will be reversed.
```bash
## Genebody
REGION="genebody"
for STRAND in both same opposite
do
    OUTPUT="${PROFDIR}${SAMPLE}.${REGION}.${STRAND}"
    ngs.plot.r -G hg38 -R $REGION -C ${MATE1REHEADER} -O $OUTPUT -P $THREADS -SS $STRAND -SE 1 -L 5000 -F chipseq -D ensembl
done

## TSS
REGION="tss"
for STRAND in both same opposite
do
    OUTPUT="${PROFDIR}${SAMPLE}.${REGION}.${STRAND}"
    ngs.plot.r -G hg38 -R $REGION -C ${MATE1REHEADER} -O $OUTPUT -P $THREADS -SS $STRAND -SE 1 -L 5000 -F chipseq -D ensembl
done
```

## TES
REGION="tes"
for STRAND in both same opposite
do
    OUTPUT="${PROFDIR}${SAMPLE}.${REGION}.${STRAND}"
    ngs.plot.r -G hg38 -R $REGION -C ${MATE1REHEADER} -O $OUTPUT -P $THREADS -SS $STRAND -SE 1 -L 5000 -F chipseq -D ensembl
done
```


#### Combining sense and anti-sense profiles
By default ngs.plot produces average profiles in .pdf format.  There is also a zip file containing the underlying data,  This may be used to combine the sense and anti-sense profiles on a single set of axes.  This is best achieved using R.


#### Scaling the profiles
A scale factor may be used to normalise for differences in library sizes across samples.  There are many ways to generate such a factor, such as the ratio of reads between two samples of spike-ins.  A more robust scale factor may be calculated using the "estimateSizeFactors" function from the Bioconductor package DESeq2 using sample gene count information.<br>

The scale factor may be used to modify the number of valid reads discovered in the BAM file given in the resulting ".cnt" file.  Edit this file to multiply read count by your scale factor and re-run ngs.plot using the same parameters.


#### Custom genes / intervals
Profiles may be restricted to specific gene sets and also to user-defined intervals.  See the ngs.plot documentation for further details.<br>
https://github.com/shenlab-sinai/ngsplot<br>
https://github.com/shenlab-sinai/ngsplot/wiki/ProgramArguments101<br>
