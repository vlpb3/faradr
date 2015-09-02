\name{CombineCountTables}
\alias{CombineCountTables}
\title{Plot number of reads in each sample for each grouping factor.}
\usage{
CombineCountTables(path, output, fileSuffix)
}
\arguments{
\item{path}{path to count tables}
\item{output}{path of the output file}
\item{fileSuffix}{suffix of the files with count tables}
}
\value{
Combined count table class (small/medium/large) of miRNA.
}
\description{
A function to read all the Index Statistics (samtools output) and combine them in one table for all samples
}

