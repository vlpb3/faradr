library(ShortRead)
library(stringr)

context("fq")

# location of sample data for testing
dataDir <- "../data"

test_that("QualityTrimReads trimmes correctly", {                            
  # test if after trimming there are no bases of quality
  # below min.qual
  min.qual <- 20
  fq.file <- file.path(dataDir, "sample.fq.gz")
  sreadq <- readFastq(fq.file) 
  
  trimmed.reads <- QualityTrimReads(sreadq, min.qual)
  
  # get all ASCII codes for qualities below min.qual
  encod <- encoding(quality(sreadq))
  q.below <- names(encod[encod %in% 1:(min.qual-1)])
  
  # fetch all qualities of trimmed reads
  quals <- quality(trimmed.reads)
  
  # iterate over all qualities to see if qual below
  # min.qual (those in q.below) are present
  high.quals.present <- FALSE
  for (i in 1:length(quals)){
    quals.string <- as.character(quals[[i]])
    if (any(str_detect(quals.string, stringr::fixed(q.below)))) {
      high.quals.present <- TRUE
      break
    }
  }
  expect_false(high.quals.present) 
})