#' Count reads in all 
#' @param bam_dir directory with alignments
#' @param base_name part of the name commot to all the alignment files
#' @param table_path file path for resulting count table
#' @importFrom stringr str_detect str_split fixed
#' @importFrom GenomicAlignments readGAlignments seqnames
#' @export
CountTable <- function(bam_dir, base_name, table_path) {
  # fetch list of bam files, with provided base name
  dir_files <- list.files(bam_dir, base_name)
  
  # ignore index files
  bam_files <- dir_files[!str_detect(dir_files, fixed(".bam.bai"))]
  
  # count alignments in each bam file
  all_counts <- lapply(bam_files, function(bam) {
    sampleid <- str_split(bam, fixed(base_name))[[1]][1]
    aln <- readGAlignments(file.path(bam_dir, bam), "BAM")
    counts <- table(as.character(seqnames(aln)))
    sample_counts <- data.frame(counts)
    names(sample_counts) <- c("seqid", sampleid)
    return(sample_counts)
  })
  
  # merge those table for each bam file inot one with samples in columns
  count_table <- Reduce(function(a, b) {
    merge(a, b, by="seqid")
  }, all_counts)
  
  # save count table into a file
  write.table(count_table, file=table_path, row.names=F,
              sep="\t", quote=F)
}
