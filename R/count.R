#' Count reads in all 
#' @param bam_dir directory with alignments
#' @param base_name part of the name commot to all the alignment files
#' @param table_path file path for resulting count table
#' @importFrom stringr str_detect
#' @export
PlotSpikeCounts <- function(bam_dir, base_name, table_path) {
  # fetch list of bam files, with provided base name
  dir_files <- list.files(bam_dir, base_name)
  
  # ignore index files
  bam_files <- dir_files[!str_detect(dir_files, fixed(".bam.bai"))]
  
  # count alignments in each bam file
  
}