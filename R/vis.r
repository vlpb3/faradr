library(ggbio)
#' Plot aligned reads over annotated region
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