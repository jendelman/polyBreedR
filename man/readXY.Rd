% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readXY.R
\name{readXY}
\alias{readXY}
\title{Read SNP array intensity data}
\usage{
readXY(filename, skip = 9, output = "ratio")
}
\arguments{
\item{filename}{filename}

\item{skip}{number of lines to skip before the header line with the column names}

\item{output}{One of three options: "ratio","theta","AD"}
}
\value{
matrix with dimensions markers x individuals
}
\description{
Read SNP array intensity data
}
\details{
The first two columns of the tab-delimited input file should be the SNP and Sample ID. Columns labeled "X" and "Y" contain the signal intensities for the two alleles. Use \code{output} to specify whether to return the ratio = Y/(X+Y) or theta = atan(Y/X)*2/pi. Option "AD" exports the XY data in the allele depth format for a VCF file ("X,Y"), with the X and Y values multiplied by 100 and rounded to the nearest integer.
}
