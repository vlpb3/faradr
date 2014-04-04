library(rtracklayer)
library(GenomicRanges)

context("vis")

# location of sample data for testing
dataDir <- "../data"
consensusSeqId <- "EMBOSS_001|consensus"

cleanup <- function() {
  files <- list.files(file.path(dataDir))
  file.pattern <- '^test_'
  test.files <- files[sapply(files,
                             function(x){str_detect(x, file.pattern)})]
  test.file.paths <- sapply(test.files,
                            function(x) {file.path(dataDir, x)})
  file.remove(test.file.paths)
}

test_that("Plotting functions in vis generates and saves a plots",{
  # input files
  gfffile <- file.path(dataDir, "rRNA.gff")
  bamfile <- file.path(dataDir, "sample.bam")
  
  # import data
  ## genome annotations
  gff <- import(gfffile, asRangedData=F)
  ## aligned reads
  gr <- GRanges(seqnames = c("EMBOSS_001|consensus"),
                ranges = IRanges(1, width=7000),
                strand = "*"
  )
  param <- ScanBamParam(which=gr)
  consensus.aln <- readGAlignments(bamfile, param=param)
  
  # get rid of unused selevels
  seqlevels(consensus.aln) <- consensusSeqId
  
  ## fetch only consensus sequence from gff
  consensus.annots <- gff[seqnames(gff)==consensusSeqId]
  plotfile <- file.path(dataDir, 'test_alnannot.png')
  
  # test PlotAlignAnnot function
  plotfile1 <- file.path(dataDir, 'test_alnannot.png')
  PlotAlignAnnot(consensus.aln, consensus.annots, plotfile1)
  expect_true(file.exists(plotfile1))
  
  # test PlotCoverageAnnot
  plotfile2 <- file.path(dataDir, 'test_covannot.png')
  PlotCoverageAnnot(consensus.aln, consensus.annots, plotfile2)  
  expect_true(file.exists(plotfile2))
})

cleanup()
