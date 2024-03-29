% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_LA.R
\name{impute_LA}
\alias{impute_LA}
\title{Impute from low to high density markers by Linkage Analysis (LA)}
\usage{
impute_LA(ped.file, high.file, low.file, low.format = "GT", out.file)
}
\arguments{
\item{ped.file}{pedigree file for progeny}

\item{high.file}{name of file with phased parental genotypes}

\item{low.file}{name of VCF file with progeny}

\item{low.format}{either "GT" (default) or "AD"}

\item{out.file}{name of CSV output file}
}
\description{
Impute from low to high density markers by Linkage Analysis
}
\details{
You must have separately installed PolyOrigin and Julia for this function to work.

The high density file contains phased parental genotypes using 0|1 format. The first 3 columns are the genetic map in cM: marker, chrom, position. To output imputed data with physical rather than genetic map positions, including a fourth column named "bp". Subsequent columns are the phased parental genotypes. 

VCF is assumed for the low-density file. The pedigree file has four columns: id, pop, mother, father, ploidy.

The output file contains the posterior maximum genotypes.

A temporary directory "tmp" is created to store intermediate files and then deleted.
}
