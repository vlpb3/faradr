\name{PlotMMFracDist}
\alias{PlotMMFracDist}
\title{Plot fraction of mismatches and indels in alignment.}
\usage{
PlotMMFracDist(bamfile, plotfile)
}
\arguments{
  \item{bamfile}{string bamfile filepath}

  \item{plotfile}{string output png filepath}
}
\description{
It is calculated by dividing read length by edit distance
(NM tag in bam file).
}

