#' Plot spike counts
#' @param counts_file csv file with spike counts
#' @export
PlotSpikes <- function(counts_file) {
  counts <- read.csv(counts_file, sep="/t")
  samples <- names(counts)[2:length(counts)]
  molten_counts <- melt(counts)
  names(molten_counts) <- c("seqid", "sample", "count")
  p <- ggplot(molten_counts, aes(seqid, count, group=sample, colour=sample))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
}

#' Plot normalised spike counts
#' @param counts_file csv file with spike counts
#' @param sample_read_couts named list, number of reads per each sample
#' @export
PlotNormSpikes <- function(counts_file, total_reads) {
  counts <- read.csv(counts_file, sep="\t")
  total_reads_df <- read.csv(total_reads, sep="\t")
  samples <- names(counts)[2:length(counts)]
  molten_counts <- melt(counts)
  names(molten_counts) <- c("seqid", "sample", "count")
  norm_counts <- molten_counts %>% group_by(sample) %>%
    transform(count_pm=count/sum(count)*0.001)
  p <- ggplot(norm_counts, aes(seqid, count_pm, group=sample, colour=sample))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
}