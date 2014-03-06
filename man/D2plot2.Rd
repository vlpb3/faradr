\name{D2plot2}
\alias{D2plot2}
\title{Plot base frequnecy for last 30 nt.}
\usage{
  D2plot2(samples)
}
\arguments{
  \item{samples}{ShortReadQ object from package ShortRead}
}
\value{
  plot object
}
\description{
  Lineplots representing frequencies of bases, one line per
  base, for last 30 nt in the read (from 3') One plot per
  sample. Samples are explected to be separage fastq files.
}

