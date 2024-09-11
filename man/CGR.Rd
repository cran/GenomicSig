\name{CGR}
\alias{CGR}
\title{
Chaos Game Representation
}
\description{
Chaos game representation is the depiction of sequence in graphical form. It converts long single dimensional sequence (in this case genetic sequence) into graphical form.
}
\usage{
CGR(data)
}
\arguments{
  \item{data}{
Input as a nucleic acid sequence of characters from fasta file.
}
}
\details{
This function produces visual image of DNA sequence different from the usual linear ordering of nucleotides.
}
\value{
This function produces a chaos game representation (CGR) plot and a data frame.
}
\references{
Deschavanne, P. J., Giron, A., Vilain, J., Fagot, G., & Fertil, B. (1999). Genomic signature: characterization and classification of species assessed by chaos game representation of sequences. Molecular biology and evolution, 16(10), 1391-1399.

Jeffrey, H. J. (1990). Chaos game representation of gene structure. Nucleic acids research, 18(8), 2163-2170.
}
\author{
Dr. Anu Sharma,
Dr. Dwijesh Chandra Mishra
}
\examples{
library(GenomicSig)
data(Genomicdata)
CGR(Genomicdata)
}