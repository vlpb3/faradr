\name{C3plot2}
\alias{C3plot2}
\title{Plot fraction of reads with partucular quality or higher per position.}
\usage{
  C3plot2(fqc, samples)
}
\arguments{
  \item{fqc}{FastQA from package ShortRead}

  \item{samples}{ShortReadQ object from package ShortRead}
}
\value{
  plot object
}
\description{
  4 line plots, for qualities 18, 20, 24 and 28, showing
  fraction of reads with one of those qualities per
  position in the read. One line per sample. Samples are
  explected to be separage fastq files.
}

