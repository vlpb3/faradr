library(stringr)

context('qc')

# location of sample data for testing
dataDir <- '../data'

cleanup <- function() {
  files <- list.files(file.path(dataDir))
  file.pattern <- '^test_'
  test.files <- files[sapply(files,
                             function(x){str_detect(x, file.pattern)})]
  test.file.paths <- sapply(test.files,
                            function(x) {file.path(dataDir, x)})
  file.remove(test.file.paths)
}

test_that("PlotMMDist outputs plot", {                           
  bamfile <- file.path(dataDir, "sample.bam") 
  plotfile <- file.path(dataDir, "test_mmplot.png")
  PlotMMDist(bamfile, plotfile)
  expect_true(file.exists(plotfile)) 
})

test_that("PlotMMFracDist outputs plot", {                           
  bamfile <- file.path(dataDir, "sample.bam") 
  plotfile <- file.path(dataDir, "test_mmpfrcplot.png")
  PlotMMFracDist(bamfile, plotfile)
  expect_true(file.exists(plotfile)) 
})

cleanup()