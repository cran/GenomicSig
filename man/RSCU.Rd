\name{RSCU}
\alias{RSCU}
\title{
Relative Synonymous Codon Usage
}
\description{
This function computes the relative frequency of each codon coding for an amino acid.
}
\usage{
RSCU(sequence)
}
\arguments{
  \item{sequence}{
Input as a nucleic acid sequence of characters from fasta file.
}
}
\details{
RSCU values are the number of times a particular codon is observed, relative to the number of times that the codon would be observed for a uniform synonymous codon usage.
}
\value{
RSCU returns the data frame with all indices.
}
\references{
Elek, A. (2018). coRdon: an R package for codon usage analysis and prediction of gene expressivity (Master's thesis, University of Zagreb. Faculty of Science. Department of Biology).
}
\author{
Dr. Anu Sharma,
Dr. Sudhir Srivastava
}
\examples{
library(GenomicSig)
RSCU(sequence= system.file("extdata/Nuc_sequence.fasta", package = "GenomicSig"))
}
