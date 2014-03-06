\name{D1plot2}
\alias{D1plot2}
\title{Plot base frequnecy for first 30 nt.}
\usage{
  D1plot2(samples)
}
\arguments{
  \item{fqc}{FastQA from package ShortRead}

  \item{samples}{ShortReadQ object from package ShortRead}
}
\value{
  plot object
}
\description{
  Linepolots representing frequencies of bases, one line
  per base, for first 30 nt in the read (from 5') One plot
  per sample. Samples are explected to be separage fastq
  files.
}

