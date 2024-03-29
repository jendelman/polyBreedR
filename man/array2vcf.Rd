% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array2vcf.R
\name{array2vcf}
\alias{array2vcf}
\title{SNP array to VCF}
\usage{
array2vcf(array.file, map.file, model.file = NULL, ploidy, vcf.file)
}
\arguments{
\item{array.file}{name of input file with SNP array allele intensities}

\item{map.file}{vcf file with map positions for the markers}

\item{model.file}{normal mixture model parameters for genotype calls}

\item{ploidy}{sample ploidy, for use with \code{model.file}}

\item{vcf.file}{output vcf file}
}
\description{
Converts output from Genome Studio (Final Report or Wide) to VCF
}
\details{
Auto-detects whether the input file is a Genome Studio Final Report, which is a "long" format with 9-row header, or in "wide" format, where all the data for each marker is one row. XY values are multiplied by 100% to generate AI field. 

Genotype calls will attempt to be imported from the GS Final Report when \code{model.file=NULL}. For diploids, columns named "Allele 1 - AB" and "Allele 2 - AB" are expected. For tetraploids, a single column named "Alleles - AB" is expected. 

It is assumed that the parameters in \code{model.file} lead to genotype calls for the dosage of allele B. For a VCF file, genotype calls need to be based on the dosage of ALT. By default, it is assumed that A is the REF allele. For variants where B is REF, include "REF=B" as INFO in the VCF \code{map.file}.
}
