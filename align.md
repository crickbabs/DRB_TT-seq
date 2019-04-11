## This bash script provides instruction for aligning paired Illumina sequence reads against a reference genome using STAR.

In this example a *S. cerevisiae* RNA spike-in was inserted into the human RNA sample prior to sequencing.  Spike-in abundance was estimated by a separate mapping of the sequence reads to the yeast genome using STAR.  However, creating a composite human/yeast genome and aligning to that or using an alignment free abundance estimator such as [kallisto](https://pachterlab.github.io/kallisto/) are also valid options.

The purpose of the spike-in is to generate a scale-factor to account for differences in library size between multiple samples.  A scale factor may be as simple as a ratio of mapped spike-in reads between two samples.  A more robust scale factor may be calculated across multiple samples using the "estimateSizeFactors" function from the Bioconductor package DESeq2 using sample specific spike-in gene count information.

---

The output of this script are two sorted and index BAM files - one for the target organism and one for the spike-in.  There are also gene level counts produced by STAR that may be may used to generate scale factors.

*It is assumed that the FASTQ reads have been checked for quality and any filtering, adapter trimming etc. that might be required has been done prior to this step.*

Dependencies:
SAMtools     http://samtools.sourceforge.net/
STAR         https://github.com/alexdobin/STAR

---

#### Set working, temporary and results directories
```bash
WORKDIR="/camp/stp/babs/working/mitterr/projects/svejstrupj/lea.gregersen/SCAF.methods_paper/work/"
TMPDIR="${WORKDIR}tmp/"
ALIGNDIR="${WORKDIR}alignments/"
SPIKEDIR="${WORKDIR}alignments_spike/"
mkdir -p $TMPDIR
mkdir -p $ALIGNDIR
mkdir -p $SPIKEDIR
```


#### Sample information: sample name, location of paired FASTQ files.
```bash
SAMPLE="WT";
FQ1="${WORKDIR}FQ1.fastq.gz"
FQ2="${WORKDIR}FQ2.fastq.gz"
#####################BAM="${ALIGNDIR}{=${SAMPLE}.bam"
```


#### Threads - set to take advantage of multi-threading and speed things up.
```bash
THREADS=8
```


#### Path to STAR genome indices
These were created using GRCh38 Ensembl v86 (*Homo sapiens*) and R64-1-1 Ensembl v86 (*Saccharomyces cerevisiae*) genome sequences and GTF files downloaded from the [Ensembl](https://www.ensembl.org/index.html) database.  Please refer to the STAR manual for information on how to create your own genomes indices.
```bash
HUMANIDX="/camp/svc/reference/Genomics/babs/homo_sapiens/ensembl/GRCh38/release-86/genome_idx/star/100bp/"
SPIKEIDX="/camp/svc/reference/Genomics/babs/saccharomyces_cerevisiae/ensembl/R64-1-1/release-86/genome_idx/rsem/star/100bp/"
```


#### Align to the human genome.  Sort and index the genome BAM.
```bash
cd $ALIGNDIR
STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${HUMANIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE} --outFileNamePrefix ${SAMPLE}.
samtools sort --threads ${THREADS} -o ${ALIGNDIR}${SAMPLE}.sorted.bam ${ALIGNDIR}${SAMPLE}.bam
samtools index ${ALIGNDIR}${SAMPLE}.sorted.bam
```


#### Align to the yeast genome (spike-in).  Sort and index the genome BAM.
```bash
cd $SPIKEDIR
STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${SPIKEIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE}.spike --outFileNamePrefix ${SAMPLE}.
samtools sort --threads ${THREADS} -o ${ALIGNDIR}${SAMPLE}.sorted.bam ${ALIGNDIR}${SAMPLE}.bam
samtools index ${ALIGNDIR}${SAMPLE}.sorted.bam
```
