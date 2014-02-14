library(ShortRead)
library(ggbio)
library(ggplot2)
library(plyr)

#' Plot aligned reads over annotated region.
#' @param aln read alignments in GAlignments objet (GenomicTRanges)
#' @param annot genome annotations in GRanges object
#' @param plotfile path for output plot
PlotAlignAnnot <- function(aln, annots, plotfile) {
    reads.track <- autoplot(aln)
    annots.track <- autoplot(annots,
                             aes(color=group, fill=group))
    png(filename=plotfile,
        width=1024, height=512)  
    print(tracks(alignments=reads.track,
                 annotation=annots.track,
                 heights=c(5,1)))  
    dev.off()
}

#' Plot aligned reads over annotated region.
#' @param aln read alignments in GAlignments objet (GenomicTRanges)
#' @param annot genome annotations in GRanges object
#' @param plotfile path for output plot
PlotCoverageAnnot <- function(aln, annots, plotfile) {
    cov <- coverage(aln)
    coverage.track <- autoplot(cov, binwidth=10)
    
    annots.track <- autoplot(annots,
                             aes(color=group, fill=group))
    png(filename=plotfile,
        width=1024, height=512)  
    print(tracks(coverage=coverage.track,
                 annotation=annots.track,
                 heights=c(5,1)))  
    dev.off()
}

#' Plot quality scores per cycle
#' @param data.dir string path to data
#' @param file.pattern string pattern for input files
#' @return plot object
PlotPerCycleQuality <- function(data.dir, file.pattern) {
    require(ggplot2)
    input.type <- "fastq"
    sreadq.qa <- qa(data.dir, file.pattern, type=input.type)
    perCycle <- sreadq.qa[["perCycle"]]
    perCycle.q <- perCycle$quality
    pcq <- perCycle.q[rep(seq_len(nrow(perCycle.q)), perCycle.q$Count), ]
    p <- ggplot(pcq, aes(factor(Cycle), Score))
    p <- p + geom_boxplot()
    return(p)
}

#' Plot cumulative qualities 
#' @param data.dir string path to data
#' @param file.pattern string pattern for input files
#' @return plot 
PlotPerCycleCumQuality <- function(data.dir, file.pattern) {
    require(ggplot2)
    require(plyr)
    fastq.files <- list.files(data.dir, file.pattern)
    q95rlen.df <- adply(fastq.files, 1, function(f, data.dir) {
      fq <- readFastq(data.dir, f)
      q95rlen <- quantile(width(fq), .95)
      data.frame(lane=c(f), q95rlen=c(q95rlen))
    }, data.dir=data.dir)
    q95rlen.df <- subset(q95rlen.df, select = c(lane, q95rlen))

    input.type <- "fastq"
    sreadq.qa <- qa(data.dir, file.pattern, type=input.type)
    perCycle <- sreadq.qa[["perCycle"]]
    pcq <- perCycle$quality
    pcq <- ddply(pcq, .(lane, Cycle), mutate, CycleCounts=sum(Count))    
    pcq$CountFrac <- pcq$Count / pcq$CycleCounts
    pcq <- ddply(pcq, .(lane, Cycle), mutate, ScoreCumSum=cumsum(CountFrac))
    pcq <- ddply(pcq, .(lane), mutate, LaneCounts=CycleCounts[1])
    pcq <- merge(x=pcq, y=q95rlen.df, by="lane", all=TRUE)
    pcq <- pcq[pcq$Cycle <= pcq$q95rlen, ]
    subpcq <- pcq[pcq$Score %in% c(15, 20, 25, 30), ]
    p <- ggplot(subpcq)
    p <- p + geom_line(aes(Cycle, 1-ScoreCumSum, group=Score, colour=factor(Score)), alpha=I(0.8))
    p <- p + geom_area(aes(Cycle, CycleCounts / LaneCounts), position="identity", alpha=I(0.1))
    p <- p + facet_wrap(~lane, scales="free_x")
    return(p)
}

#' Plot mean quality per read
#' @param data.dir path to data dir
#' @param file.pattern pattern of input fastq
#' @return plot object
PlotMeanReadQuality <- function(data.dir, file.pattern) {
  count.qmeans <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    qm <- as(quality(fq), "matrix")
    row.means <- rowMeans(qm, na.rm=T)
    qmeans <- data.frame(mean.qual=row.means, lane=rep(s, length(row.means)))
    return(qmeans)
  }
  samples <- list.files(data.dir, file.pattern)
  qmeans <- adply(samples, 1, count.qmeans, data.dir=data.dir)
  p <- ggplot(qmeans, aes(x=mean.qual))
  p <- p + geom_histogram(aes(y=..density..),
                          alpha=.1,
                          fill="green",
                          colour="darkgreen",
                          binwidth=.5)
  p <- p + geom_density() + facet_wrap(~lane)
  return(p)
}

#' Plot read length distribution
#' @param data.dir path to data dir
#' @param file.pattern pattern of input fastq
#' @return plot object 
PlotReadLengthDistribution <- function(data.dir, file.pattern) {
  samples <- list.files(data.dir, file.pattern)
  count.rlens <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    rlens <- width(fq)
    rlen95 <- quantile(width(fq), .95)
    rlens <- rlens[rlens <= rlen95]
    rlens.df <- data.frame(read.len=rlens, lane=rep(s, length(rlens)))
    return(rlens.df)
  }
  
  rlens <- adply(samples, 1, count.rlens, data.dir=data.dir)
  p <- ggplot(rlens, aes(x=read.len))
  p <- p + geom_histogram(aes(y=..density..),
                          alpha=0.02,
                          fill="green",
                          colour="darkgreen",
                          binwidth=1)
  p <- p + geom_density(adjust=3) + facet_wrap(~lane)
  return(p)
}

#' Plot base call per cycle
#' @param data.dir string path to data
#' @param file.pattern string pattern for input files
#' @param type stirng input type, 'fastq' is default value
#' @return plot 
PlotPerCycleBaseCalls <- function(data.dir, file.pattern, type="fastq") {
    require(ggplot2)  
    require(plyr)
    
    sreadq.qa <- qa(data.dir, file.pattern, type=type)
    perCycleBaseCall <- sreadq.qa[["perCycle"]]$baseCall
    perCycleCounts  <- ddply(perCycleBaseCall, c("Cycle"), summarise,
                           countsInCycle = sum(Count))
    perCycleBaseCall$callsFraction <- perCycleBaseCall$Count / perCycleCounts$countsInCycle[perCycleBaseCall$Cycle] 
    p <- ggplot(perCycleBaseCall, aes(Cycle, callsFraction))
    p <- p + geom_point(aes(colour=factor(Base)), alpha=0.7)
    p <- p + geom_line(aes(colour=factor(Base)), alpha=0.7) 
    return(p)
}