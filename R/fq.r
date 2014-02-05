library(assertthat)
library(ShortRead)
#' Trim fastq reads.
#' 
#' Trim reads from 3' end leaving only bases of quality
#' equal or above min.quality.
#' @param sreadq ShortReadQ object
#' @param min.quality integer minimum allowed quality
#' @return ShortReadQ object with trimmed reads
QualityTrimReads <- function(sreadq, min.quality) {
    
    # guess encoding from the reads
    encod <- encoding(quality(sreadq))
    ascii.quality <- names(encod[encod == (min.quality-1)])
    trimmed.sreadq <- trimTails(sreadq, 1, ascii.quality, successive=F)

    return (trimmed.sreadq)
}


