library(ShortRead)
library(ggbio)
library(ggplot2)
library(plyr)

#' Plot aligned reads over annotated region.
#' @param aln read alignments in GAlignments objet (GenomicTRanges)
#' @param annot genome annotations in GRanges object
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
#' @param type stirng input type, 'fastq' is default value
#' @return plot 
PlotPerCycleQuality <- function(data.dir, file.pattern, type="fastq") {
    require(ggplot2)
    max.rlen <- 200
    sreadq.qa <- qa(data.dir, file.pattern, type=type)
    perCycle <- sreadq.qa[["perCycle"]]
    perCycle.q <- perCycle$quality
    perCycle.q <- perCycle.q[perCycle.q$Cycle <= max.rlen, ]
    pcq <- perCycle.q[rep(seq_len(nrow(perCycle.q)), perCycle.q$Count), ]
    p <- ggplot(pcq, aes(factor(Cycle), Score))
    p <- p + geom_boxplot()
    return(p)
}

#' Plot mean quality per read
#' @param sreadq ShortReadQ object
#' @return plot
PlotMeanReadQality <- function(sreadq) {
    qm <- as(quality(sreadq), "matrix") 
    qmeans <- data.frame(mean_quality=rowMeans(qm, na.rm=T))
    p <- qplot(mean_quality, data=qmeans, geom="histogram",
               binwidth=0.5, alpha=I(0.7)) 
    return(p)
}

#' Plot read length distribution
#' @param sreadq ShortReadQ object
#' @return plot
PlotReadLengthDistribution <- function(sreadq) {
    max.rlen <- 75
    rlens <- data.frame(read_len=width(sreadq[width(sreadq) <= 75]))
    p <- qplot(read_len, data=rlens, geom="histogram",
               binwidth=1, alpha=I(0.7))
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

