library(GenomicRanges)
library(ShortRead)
library(Rsamtools)

#' Plot distribution of mismatches and indels in alignment. 
#' 
#' @param bamfile string bamfile filepath
#' @param plotfile string output png filepath 
PlotMMDist <- function(bamfile, plotfile) {
  param <- ScanBamParam(tag="NM")
  aln <- readGAlignments(bamfile, param=param)
  mm.counts <- table(mcols(aln)$NM)
  png(filename=plotfile, width=1024, height=512)
  barplot(mm.counts,
          main="Mismatch distribution in the alignments",
          xlab="Number of mismatches and indels within a read alignment",
          ylab="Number of reads",
          cex.lab=2, cex.main=2)
  dev.off()
}

#' Plot fraction of mismatches and indels in alignment.
#'
#' It is calculated by dividing read length
#' by edit distance (NM tag in bam file).
#' @param bamfile string bamfile filepath
#' @param plotfile string output png filepath 
PlotMMFracDist <- function(bamfile, plotfile) {
  param <- ScanBamParam(tag="NM")
  aln <- readGAlignments(bamfile, param=param)
  mcols(aln)$mmfrac <- round(mcols(aln)$NM / width(aln), 2)
  breaks <- max(mcols(aln)$mmfrac) * 100
  png(filename=plotfile, width=1024, height=512)
  hist(mcols(aln)$mmfrac,
       breaks=breaks,
       main="Mismatch distribution in the alignments",
       xlab="Number of mismatches and indels within a read alignment",
       ylab="Number of reads",
       cex.lab=2, cex.main=2)
  dev.off()
}

#' Get quality statistics.
#' 
#' Calculates number of bases ans average lenght of the read
#' with quality above provided.
#' @param fastq stiring path to fastq file
#' @param q.threshold integer minimum allowed quality
#' @return list with average read length and gnumber of bases
#'  meeting quality requirements
GeneralQualityStats <- function(fastq, q.threshold) {
  encod <- encoding(quality(fastq)) 
  if (! q.threshold %in% encod) {
    min.qual <- range(encod)[1]
    max.qual <- range(encod)[2]
    print(paste("Quality threshold has to be within the range from",
                min.qual, "to", max.qual))   
  }
  q <- names(encod[encod == q.threshold])
  trimmed.fq <- trimTails(fastq, 1, q, successive=F)  
  mean.rlen <- mean(width(trimmed.fq)) 
  n.qual.bases <- sum(width(trimmed.fq))  
  return(list(mean.rlen=mean.rlen, n.qual.bases=n.qual.bases))
}

#' Make read statistics table
#' 