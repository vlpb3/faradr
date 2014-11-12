#' Plot spike counts
#' @param counts_file csv file with spike counts
#' @importFrom ggplot2 ggplot geom_point geom_line theme
#' @export
PlotSpikeCounts <- function(counts_file) {
  counts <- read.csv(counts_file, sep="\t")
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
#' @importFrom ggplot2 ggplot geom_point geom_line theme
#' @importFrom reshape2 melt
#' @importFrom dplyr %>% group_by
#' @export
PlotNormalSpikeCounts <- function(counts_file, total_reads) {
  # load count tables into dataframes
  counts <- read.csv(counts_file, sep="\t")
  colnames <- c("sample", "total")
  total_reads_df <- read.csv(total_reads, header=F, col.names=colnames, sep="\t")
  
  samples <- names(counts)[2:length(counts)]  

  molten_counts <- melt(counts)
  names(molten_counts) <- c("seqid", "sample", "count")
  molten_counts$total <- rep(0, nrow(molten_counts))
  for (s in samples){
    total <- total_reads_df[total_reads_df$sample == s,]$total
    molten_counts[molten_counts$sample == s,]$total <- total
  }
    
  # normalize per milion of reads (total)
  norm_counts <- molten_counts %>% group_by(sample) %>%
    transform(count_pm=(count*1e6)/total)

  p <- ggplot(norm_counts, aes(seqid, count_pm, group=sample, colour=sample))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
}

#' Normalise counts based on normalisation spike counts
#' @param spike_counts file with count table of normalization spikes
#' @param counts file with counts table to normalize
#' @param norm_counts file that will be used for calculated normalized counts
#' @export
normalize_counts <- function(spike_counts, counts, norm_counts) {
  get_norm_factors <- function(counts, locfunc=median)  {
    loggeomeans <- rowMeans(log(counts))
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans)]))
    })
  }
  spike_df <- read.table(spike_counts, row.names=1)
  norm_factors <- get_norm_factors(spike_df)
  counts_df <- read.table(counts, row.names=1, header = T)
  norm_df <- t(apply(counts_df,1, function(r){r / norm_factors}))
  write.table(norm_df, norm_counts, sep = "\t", quote=F)
}
