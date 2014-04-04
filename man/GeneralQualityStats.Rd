\name{GeneralQualityStats}
\alias{GeneralQualityStats}
\title{Get quality statistics.}
\usage{
GeneralQualityStats(fastq, q.threshold)
}
\arguments{
  \item{fastq}{stiring path to fastq file}

  \item{q.threshold}{integer minimum allowed quality}
}
\value{
list with average read length and gnumber of bases meeting
quality requirements
}
\description{
Calculates number of bases ans average lenght of the read
with quality above provided.
}

