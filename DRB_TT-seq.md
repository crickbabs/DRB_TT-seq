This is a companion script to the publication below.  It describes a pipeline for calling RNA Pol II transcription wave peak positions and elongation rates from DRB/TT-seq time-series data using R.  Instructions are given for calculating wave peaks at both the single-gene and meta-gene level.  
  

*Nascent transcriptome profiles and measurement of transcription elongation using TT-seq.*
Lea H. Gregersen<sup>1</sup> Richard Mitter<sup>2</sup> and Jesper Q. Svejstrup<sup>1</sup>
<sup>1</sup>Mechanisms of Transcription Laboratory, The Francis Crick Institute, 1 Midland Road, London, NW1 1AT, UK
<sup>2</sup>Bioinformatics and Biostatistics, The Francis Crick Institute, 1 Midland Road, London NW1 1AT, UK

---


#### Define genic intervals and calculate coverage directly from BAM files.

```bash
library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(bamsignals)
library(reshape2)
library(ggplot2)
library(DT)
```



```bash
work.dir <- "/path/to/my/working_directory/"
setwd(work.dir)
```


#### Create a data.frame detailing BAM files and sample names and any pertinent meta-data such as DRB release time.
```bash
drb.df <- data.frame(
  name = c("DRB.10min","DRB.20min","DRB.30min","DRB.40min"),
  bam  = c("bam/DRB.10min.bam","bam/DRB.20min.bam","bam/DRB.30min.bam","bam/DRB.40min.bam"),
  time = c(10,20,30,40),
  stringsAsFactors=F)
drb.df <- drb.df[order(drb.df$time),]
```


#### Genome transcriptome GTF files are available for download from [Ensembl](http://www.ensembl.org/index.html).  The GTF file may be indexed using [igvtools](https://software.broadinstitute.org/software/igv/igvtools_commandline) *index* command. Make sure that the assembly matches the one used for read alignments.
```bash
gtf.file     <- "data/Homo_sapiens.GRCh38.86.gtf"
```


#### Filters genes to remove any that are overly short, overly long, any non-coding genes any overlapping another gene and any that map to non-standard chromosomes.
```bash
gtf.dat <- import(gtf.file)
gene.gr <- gtf.dat[gtf.dat$type %in% "gene",]
names(gene.gr) <- gene.gr$gene_id
gene.gr$gene.width <- width(gene.gr)

gene.gr <- gene.gr[width(gene.gr) >= 60000 & width(gene.gr) < 300000]
gene.gr <- gene.gr[gene.gr$gene_biotype %in% "protein_coding"]
gene.gr <- gene.gr[countOverlaps(gene.gr,gene.gr) == 1]
gene.gr <- keepSeqlevels(gene.gr,c(as.character(1:22),"X","Y"),pruning.mode="coarse")
```


#### Create GRanges object representing filtered genes -2kb:+120kb around their promoters.
```bash
up.ext  <- 2000
dn.ext  <- 120000
intervals.gr <- promoters(gene.gr, upstream=up.ext, downstream=dn.ext)
intervals.gr <- trim(intervals.gr)
intervals.gr <- intervals.gr[width(intervals.gr) == (up.ext+dn.ext)]
```


#### Calculate bp-level coverage over the extended gene intervals.
Currently this is not strand-specific as "bamCoverage" doesn't allow it, but filtering out overlapping genes should take care of most of the problems.  Alternatively, one could filter the bam file by strand prior to calculating coverage.
```bash
sigs.list <- list()
for (r in 1:nrow(drb.df)) {
    sigs <- bamCoverage(
        bampath                    = drb.df$bam[r],
        gr                         = intervals.gr,
        filteredFlag               = 1024,  # remove duplicates
        paired.end                 = "ignore",
        mapq                       = 20,
        verbose                    = FALSE)
    sigs.list[[drb.df$name[r]]] <- t(alignSignals(sigs))
    rownames(sigs.list[[drb.df$name[r]]]) <- names(intervals.gr)
}
```


#### Scale to Read counts Per Million (RPM).
```bash
read.length <- 75
rpm.list   <- list()
for (n in 1:length(sigs.list)) {
    n.sig <- sigs.list[[n]] / read.length
    n.sum <- sum(n.sig)
    n.sf  <- n.sum / 1000000
    rpm.list[[names(sigs.list)[n]]] <- n.sig / n.sf
}
```

***

## Meta-gene level wave peak calling.

#### Create meta-profiles by taking a trimmed mean.
```bash
meta.dat <- sapply(rpm.list,function(x){ apply(x,2,function(y) { mean(y,na.rm=T,trim=0.01) }) })
```


#### Fit a smoothing spline to each meta-profile.
```bash
# The spar parameter might need tweaking depending on the fit
spline.dat <- t(apply(meta.dat,2,function(x){smooth.spline(1:length(x),x,spar=0.9)$y }))
```


#### Plot the meta-profiles.
```bash
# Coverage
cov.plot <- melt(meta.dat)
colnames(cov.plot) <- c("position","name","RPM")
cov.plot$position  <- cov.plot$position-up.ext-1
cov.plot$time      <- drb.df$time[match(cov.plot$name,drb.df$name)]

# Spline
spline.plot <- melt(t(spline.dat))
colnames(spline.plot) <- c("position","name","RPM")
spline.plot$position  <- spline.plot$position-up.ext-1
spline.plot$time      <- drb.df$time[match(spline.plot$name,drb.df$name)]

P1 <- ggplot(cov.plot,aes(x=position,y=RPM,colour=name)) + geom_line(alpha=0.4)
P1 <- P1 + geom_line(aes(x=position,y=RPM,colour=name),spline.plot,size=2)
P1 <- P1 + xlab("Position relative to TSS (kb)") + ylab("RPM") + ggtitle("DRB/TT-seq coverage")
P1 <- P1 + scale_x_continuous(breaks=c(0,40000,80000,120000),label=c("TSS","40kb","80kb","120kb"))
#P1 <- P1 + scale_colour_manual(values=c("DRB.10min"="#231F20", "DRB.20min"="#58595B", "DRB.30min"="#A7A9AC", "DRB.40min"="#D1D3D4"))
#ggsave(filename="results/DRB_metaprofile.png",plot=P1,device="png",height=5)
P1
```


#### Calculate wave peaks from meta-profiles as the maximum point on the spline.
Maxima are only called after the preceeding timepoint's maximum position, forcing the wave to advance with time.
```bash
wf.dat <- drb.df[,c("name","time")]
wf.dat$wave <- 0
for (r in 1:nrow(wf.dat)) {
  present.sample <- wf.dat$name[r]
  if (r==1) {
    previous.wf <- 0
    wf.dat$wave[r] <- which.max(spline.dat[present.sample,])
  } else {
    previous.wf    <- wf.dat$wave[r-1]
    wf.dat$wave[r] <- which.max(spline.dat[present.sample,-1:-previous.wf])+previous.wf
  }
}
wf.dat$wave <- wf.dat$wave - 2001
wf.dat$wave <- wf.dat$wave / 1000

# Optionally add an additional time=0 datapoint which assumes a wave peak at position=0.
wf.dat <- data.frame(rbind(c("DRB.0min",0,0),wf.dat),stringsAsFactors=F)
wf.dat$time <- as.numeric(wf.dat$time)
wf.dat$wave <- as.numeric(wf.dat$wave)
```


#### Wave peaks calculated from the meta-gene profiles.
```bash
datatable(wf.dat,rownames=FALSE)
```


#### Fit a linear model to the wave peak positions as a function of time to determine the rate of elongation, kb/min.
```bash
lm.fit <- lm(wf.dat$wave~wf.dat$time)
elongation.rate <- lm.fit$coefficients[2]
P2 <- ggplot(wf.dat,aes(x=time,y=wave)) + geom_point(size=4)
P2 <- P2 + xlab("Time (min)") + ylab("Wave position (kb)")
P2 <- P2 + geom_abline(intercept = lm.fit$coefficients[1], slope = lm.fit$coefficients[2])
P2 <- P2 + annotate(geom="text", x=10, y=75, label=paste("y = ",round(lm.fit$coefficients[2],2),"x",round(lm.fit$coefficients[1],2),sep=''))
#ggsave(filename="results/DRB_metaprofile_elongation_rate.png",plot=P2,device="png")
P2
```


***

## Single gene level wave peak calling

#### Generate a set of gene ids that pass an arbitrary expression threshold.
```bash
expr.mat  <- sapply(rpm.list,function(x){apply(x,1,function(y){sum(y,na.rm=T)})})
expr.plot <- ggplot(melt(expr.mat),aes(x=log2(value+0.1),fill=Var2)) + geom_histogram(bins=50)
expr.plot <- expr.plot + geom_vline(aes(xintercept=log2(100+1)),colour="darkred") + facet_grid(Var2~.) + xlab("log2(RPM+0.1)") + ylab("frequency") + ggtitle("Expression filter")
expr.gids <- rownames(expr.mat)[rowSums(expr.mat > 100) == ncol(expr.mat)]
#ggsave(filename="results/DRB_expression_filter.png",plot=expr.plot ,device="png")
expr.plot
```


#### Fit a smoothing spline over the RPM data for each gene for each sample.
```bash
spline.list <- list()
for (n in 1:length(rpm.list)) {
    spline.dat <- t(apply(rpm.list[[n]],1,function(x){smooth.spline(1:length(x),x,spar=0.9)$y }))
    spline.list[[names(rpm.list)[n]]] <- spline.dat
}
```


#### Calculate wave peak for each gene as the maximum point on the spline.
```bash
wf.genes <- sapply(spline.list,function(x){ apply(x,1,function(y){which.max(y)}) })
rownames(wf.genes) <- rownames(spline.list[[1]])
```


#### Filter the gene level wave peak predictions.  Remove any genes that are lowly expressed, have missng values, have duplicate values or whose wave doesn't advance with time.  Select only genes with a wave-peak after the first 2 kb in the 10min sample.
```bash
wf.genes.filt <- wf.genes[expr.gids,c("DRB.10min","DRB.20min","DRB.30min","DRB.40min")]
wf.genes.filt <- wf.genes.filt[!apply(wf.genes.filt,1,function(x){any(is.na(x))}),]
wf.genes.filt <- wf.genes.filt[!apply(wf.genes.filt,1,function(x){any(duplicated(x))}),]
wf.genes.filt <- wf.genes.filt[apply(wf.genes.filt,1,function(x){all(order(x)==1:nrow(drb.df))}),]
wf.genes.filt <- wf.genes.filt[wf.genes.filt[,"DRB.10min"]>2000,]

```


#### Sometimes it is necessary to disregard the final timepoint when generating the filter as transcription might have already reached the end of the gene at that point.  This version only uses the first three timepoints.
```bash
# wf.genes.filt <- wf.genes[expr.gids,c("DRB.10min","DRB.20min","DRB.30min")]
# wf.genes.filt <- wf.genes.filt[!apply(wf.genes.filt,1,function(x){any(is.na(x))}),]
# wf.genes.filt <- wf.genes.filt[!apply(wf.genes.filt,1,function(x){any(duplicated(x))}),]
# wf.genes.filt <- wf.genes.filt[apply(wf.genes.filt,1,function(x){all(order(x)==1:3)}),]
# wf.genes.filt <- wf.genes.filt[wf.genes.filt[,"DRB.10min"]>2000,]

```


#### Compile results.
```bash
wf.genes.dBase <- data.frame(
    gene_id     = rownames(wf.genes),
    gene_name   = intervals.gr[rownames(wf.genes),]$gene_id,
    gene_width  = intervals.gr[rownames(wf.genes),]$gene.width,
    mean.rpkm   = rowMeans(expr.mat[rownames(wf.genes),]),
    wf.order    = apply(wf.genes[rownames(wf.genes),drb.df$name],1,function(x){paste(order(x),sep='',collapse='')}),
    wf.genes,
    filter      = rownames(wf.genes) %in% rownames(wf.genes.filt),
    stringsAsFactors=F,
    check.names=F)
```


#### Fit a linear model to the wave peak positions as a function of time to determine the rate of elongation, kb/min.
```bash
lm.dat <- wf.genes.dBase[,drb.df$name]
wf.genes.dBase$elongation.rate <- apply(lm.dat,1,function(x) {
  time <- c(0,drb.df$time)
  wf   <- c(0,x)/1000
  lm(wf~time)$coef[2]
})
```


#### Plot a distribution of elongation rates for genes passing the filter.
```bash
single.elong_rate.dat <- data.frame(wf.genes.dBase[wf.genes.dBase$filter,c("gene_name","elongation.rate")],stringsAsFactors=F)
P3 <- ggplot(single.elong_rate.dat,aes(x=elongation.rate)) + geom_histogram(alpha=0.5,colour="#7CAE00",fill="#7CAE00")
P3 <- P3 + ylab("frequency") + xlab("Elongation rate (kb/min)") + ggtitle(paste("Elongation rate, n=",nrow(single.elong_rate.dat),sep=''))
#ggsave(filename="results/DRB_elongation_rate_distribution.png",plot=P3,device="png")
P3
```


#### Session information
```bash
sessionInfo()
```
