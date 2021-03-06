% Generated by roxygen2 (4.0.2): do not edit by hand
\name{C3.design.plot}
\alias{C3.design.plot}
\title{Plot fraction of reads with partucular quality or higher per position.}
\usage{
C3.design.plot(samples, fqc, design.table)
}
\arguments{
\item{samples}{ShortReadQ object from package ShortRead}

\item{fqc}{FastQA from package ShortRead}

\item{design.table}{data.frame holds information about experimantal design}
}
\value{
list of plot objects
}
\description{
4 line plots, for qualities 18, 20, 24 and 28, showing fraction
of reads with one of those qualities per position in the read.
One line per sample.
Samples are explected to be separage fastq files.
There is one plot per grouping factor.
Samples colored by grouping factor.
Grouping factors are extracted from design.table.
}

