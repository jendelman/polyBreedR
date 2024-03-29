% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dart2vcf.R
\name{dart2vcf}
\alias{dart2vcf}
\title{Convert DArTag to VCF}
\usage{
dart2vcf(counts.file, dosage.file, vcf.file, ploidy, first.data.row = 9)
}
\arguments{
\item{counts.file}{DArTag collapsed counts file}

\item{dosage.file}{DArTag dosage file}

\item{vcf.file}{name of VCF output file}

\item{ploidy}{ploidy}

\item{first.data.row}{default is 9 for DArTag format}
}
\description{
Convert DArTag to VCF
}
\details{
Two input files expected. counts.file is the two-row collapsed counts file, whereas dosage.file has one row per target, with chrom and position in columns 4 and 5. DArT reports dosage of REF, whereas VCF standard is based on dosage of ALT. The dosage is exported as GT field in VCF.

Duplicate samples are renamed by appending the "Target ID".
}
