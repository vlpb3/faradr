library(ShortRead)
library(ggbio)
library(ggplot2)

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

#' Plot qaality scores per cycle
#' @param data.dir string path to data
#' @param file.pattern string pattern for input files
#' @param type stirng input type, 'fastq' is default value
#' @return plot 
PlotPerCycleQualityScores <- function(data.dir, file.pattern, type="fastq") {
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
