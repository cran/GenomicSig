\name{GC_content}
\alias{GC_content}
\title{
GC content of nucleic acid sequences
}
\description{
This function calculates the percentage of guanine (G) or cytosine (C) nitrogenous bases in a DNA or RNA molecule. This measure indicates the proportion of G and C bases out of an implied four total bases, which including adenine (A) and thymine (T) in DNA. And adenine and uracil (U) in RNA along with G and C.
}
\usage{
GC_content(sequence)
}
\arguments{
  \item{sequence}{
Input as a nucleic acid sequence of characters from fasta file.
}
}
\details{
G+C content is estimated with ambiguous bases taken into account.
}
\value{
This function returns the fraction of G+C as a numeric vector of length one for all sequences.
}
\references{
Nakabachi, A., Yamashita, A., Toh, H., Ishikawa, H., Dunbar, H. E., Moran, N. A., & Hattori, M. (2006). The 160-kilobase genome of the bacterial endosymbiont Carsonella. Science, 314(5797), 267-267.

Barbu, E., Lee, K. Y., & Wahl, R. (1956, August). Content of purine and pyrimidine base in desoxyribonucleic acid of bacteria. In Annales de l'Institut Pasteur (Vol. 91, No. 2, p. 212).
}
\author{
Dr. Anu Sharma,
Dr. Girish Kumar Jha
}
\examples{
library(GenomicSig)
GC_content(sequence = system.file("extdata/Nuc_sequence.fasta", package = "GenomicSig"))
}
