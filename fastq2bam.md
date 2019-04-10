BLAH

```bash
#!/usr/bin/bash 

WORKDIR="/camp/stp/babs/working/mitterr/projects/svejstrupj/lea.gregersen/SCAF.methods_paper/work/"
TMPDIR="${WORKDIR}tmp/"
mkdir -p $TMPDIR
mkdir -p $BIGWIGDIR
```

It might be necessary to re-header the BAM file so that the chromosome names match those in the ngs.plot database, e.g. standard chromosomes preceeded with "chr" for hg38.
This is most easily achieved using SAMtools' "reheader" function.
For human hg38 Ensembl alignments the following steps should do the job.
This step may be skipped if your header is already compliant.
```bash
MATE1REHEADER="${WORKDIR}WT.mate1.reheader.bam"
samtools view -H ${MATE1} | sed -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - ${MATE1} > ${MATE1REHEADER}
samtools index ${MATE1REHEADER}
```

