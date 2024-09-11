\name{AMIP}
\alias{AMIP}
\title{
Average Mutual Information Profile (AMIP)
}
\description{
The Average Mutual Information Profile (AMIP) detects long-range correlations in a given DNA sequence by estimating the shared information between nucleotides situated k bases apart.
}
\usage{
AMIP(fasta_file, n1 = 1, n2 = 4)
}
\arguments{
  \item{fasta_file}{
Path to the input FASTA file containing the DNA sequence.
}
  \item{n1}{
The starting position (in bases) for Mutual Information calculation.
}
  \item{n2}{
The end position (in bases) for Mutual Information calculation.
}
}
\details{
The Average Mutual Information (AMI) provides a statistical estimate of the shared information between nucleotides situated k bases apart in the DNA sequence, where k ranges from \code{n1} to \code{n2}. This method helps identify potential patterns or correlations in the nucleotide arrangement.
}
\value{
This function returns a data frame containing the mutual information values for the specified nucleotide positions.
}
\references{
Bauer, M., Schuster, S. M., & Sayood, K. (2008). The average mutual information profile as a genomic signature. BMC Bioinformatics, 9(1), 48.
}
\author{
Dr. Anu Sharma,
Dr. Shashi Bhushan Lal
}
\examples{
library(GenomicSig)
AMIP(fasta_file = system.file("extdata/Nuc_sequence.fasta", package = "GenomicSig"), n1 = 1, n2 = 4)
}
