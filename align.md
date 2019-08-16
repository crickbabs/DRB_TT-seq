## This bash script provides instruction for aligning paired Illumina sequence reads against a reference genome using STAR.

In this example a *S. cerevisiae* RNA spike-in was inserted into the human RNA sample prior to sequencing.  Spike-in abundance was estimated by a separate mapping of the sequence reads to the yeast genome using STAR.  However, creating a composite human/yeast genome and aligning to that or using an alignment free abundance estimator such as [kallisto](https://pachterlab.github.io/kallisto/) are also valid options.

The purpose of the spike-in is to generate a scale-factor to account for differences in library size between multiple samples.  A scale factor may be as simple as a ratio of mapped spike-in reads between two samples.  A more robust scale factor may be calculated across multiple samples using the "estimateSizeFactors" function from the [Bioconductor](https://bioconductor.org/) package [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) using sample specific spike-in gene count information.

---

The output of this script are two sorted, duplicate marked and indexed BAM files - one for the target organism and one for the spike-in.  There are also gene level counts produced by STAR ("*.ReadsPerGene.out.tab") that may be may used to generate scale factors.

*It is assumed that the FASTQ reads have been checked for quality and any filtering, adapter trimming etc. that might be required has been done prior to running this script.*

If there are multiple sets of FASTQ files per sample, i.e. more than 1 set of paired reads, it is recommended to align these separately and merge the resultant BAM files using "samtools merge" before continuing with the downstream analysis.

Marking duplicate reads isn't strictly necessary but the calculated levels of duplication provide insight into sample quality.  Also, counts of unique mapped reads may be useful for scale factor normalisation.

Dependencies:<br>
SAMtools     http://samtools.sourceforge.net/<br>
Picard       https://broadinstitute.github.io/picard/<br>
STAR         https://github.com/alexdobin/STAR<br>

---

#### Set working, temporary and results directories
```bash
WORKDIR="/path/to/my/working_directory/"
TMPDIR="${WORKDIR}tmp/"
ALIGNDIR="${WORKDIR}alignments/"
SPIKEDIR="${WORKDIR}alignments_spike/"
mkdir -p $TMPDIR
mkdir -p $ALIGNDIR
mkdir -p $SPIKEDIR
```


#### Sample information: sample name and location of paired FASTQ files.
```bash
SAMPLE="WT";
FQ1="${WORKDIR}FQ1.fastq.gz"
FQ2="${WORKDIR}FQ2.fastq.gz"
```


#### Threads - set to take advantage of multi-threading and speed things up.
```bash
THREADS=8
```


#### Path to STAR genome indices
These were created using GRCh38 Ensembl v86 (*Homo sapiens*) and R64-1-1 Ensembl v86 (*Saccharomyces cerevisiae*) genome sequences and GTF files downloaded from the [Ensembl](https://www.ensembl.org/index.html) database.  Please refer to the STAR manual for information on how to create your own genomes indices.
```bash
HUMANIDX="/path/to/my/human_genome_index/"
SPIKEIDX="/path/to/my/yeast_genome_index/"
```


#### Align to the human genome.  Sort, mark duplicates and index the genome BAM.
```bash
cd $ALIGNDIR
STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${HUMANIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE} --outFileNamePrefix ${ALIGNDIR}${SAMPLE}.
samtools sort --threads ${THREADS} -o ${ALIGNDIR}${SAMPLE}.sorted.bam ${ALIGNDIR}${SAMPLE}.Aligned.out.bam
java -jar picard.jar MarkDuplicates INPUT=${ALIGNDIR}${SAMPLE}.sorted.bam OUTPUT=${ALIGNDIR}${SAMPLE}.sorted.marked.bam METRICS_FILE=${ALIGNDIR}${SAMPLE}.sorted.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=${TMPDIR}${SAMPLE}
samtools index ${ALIGNDIR}${SAMPLE}.sorted.marked.bam
rm ${ALIGNDIR}${SAMPLE}.sorted.bam
```

#### Align to the yeast genome (spike-in).  Sort, mark duplicates and index the genome BAM.
```bash
cd $SPIKEDIR
STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${SPIKEIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE}.spike --outFileNamePrefix ${ALIGNDIR}${SAMPLE}.
samtools sort --threads ${THREADS} -o ${SPIKEDIR}${SAMPLE}.sorted.bam ${SPIKEDIR}${SAMPLE}.Aligned.out.bam
java -jar picard.jar MarkDuplicates INPUT=${SPIKEDIR}${SAMPLE}.sorted.bam OUTPUT=${SPIKEDIR}${SAMPLE}.sorted.marked.bam METRICS_FILE=${SPIKEDIR}${SAMPLE}.sorted.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=${TMPDIR}${SAMPLE}
samtools index ${SPIKEDIR}${SAMPLE}.sorted.marked.bam 
rm ${SPIKEDIR}${SAMPLE}.sorted.bam
```
