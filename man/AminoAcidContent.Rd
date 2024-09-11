\name{AminoAcidContent}
\alias{AminoAcidContent}
\title{
Amino Acid Content
}
\description{
Amino acid content refers to the relative frequencies of amino acids used in a protein or a proteome with 20 different amino acids as 20 dimensional vectors.
}
\usage{
AminoAcidContent(fasta_file, type = c("DNA", "protein"))
}
\arguments{
  \item{fasta_file}{
Path of a fasta file containing nucleotide or protein sequence.
}
\item{type}{
Type of the sequence, can be either "DNA" or "protein".
}
}
\details{
Amino acid content refers to the relative frequencies of amino acids in the protein.If DNA sequence is given as input, it will be translated to a protein sequence. Then amino acid content will be calculated.
}
\value{
This function returns a data frame containing the sequence identifier fetched from the input fasta file and amino acid content of that sequence.
}
\references{
Sandberg, R., Bränden, C. I., Ernberg, I., & Cöster, J. (2003). Quantifying the species-specificity in genomic signatures, synonymous codon choice, amino acid usage and G+ C content. Gene, 311, 35-42.

Xiao, N., Cao, D. S., Zhu, M. F., & Xu, Q. S. (2015). protr/ProtrWeb: R package and web server for generating various numerical representation schemes of protein sequences. Bioinformatics, 31(11), 1857-1859.
}
\author{
Dr. Anu Sharma,
Megha Garg
}
\examples{
library(GenomicSig)
AminoAcidContent(fasta_file= system.file("extdata/Nuc_sequence.fasta", package = "GenomicSig"),
type = "DNA")
AminoAcidContent(fasta_file= system.file("extdata/prot_sequence.fasta", package = "GenomicSig"),
type = "protein")
}
