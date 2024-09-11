\name{OFDEG}
\alias{OFDEG}
\title{
Oligonucleotide Frequency Derived Error Gradient
}
\description{
Oligonucleotide Frequency Derived Error Gradient computes approximate convergence rate of oligonucleotide frequencies with subsequent increasing sequence length.
}
\usage{
OFDEG(sequence, c, rc, d, m, t, k, norm=0)
}
\arguments{
  \item{sequence}{
Input is a fasta file nucleic acid sequence.It accepts RData object of the fasta file
}
  \item{c}{
Minimum sequence cutoff c (which corresponds to the length of the shortest sequence in the
data set). Default is 160.
}
  \item{rc}{
Cutoff of Resampling Depth (Number of subsequence of cutoff length). Default is set to 10.
}
  \item{d}{
Sampling depth (The sampling depth refers to the number of equal length sub-sequences randomly selected from the
entire sequence). Default is set to 10. Larger sequence lengths will require greater sampling depths.
}
  \item{m}{
Word size which is initial subsequence length. Default is set to 100.
}
  \item{t}{
Step size (The step size is the change in sub-sequence length from one sampling instance to the next). Default is set to 6.
}
  \item{k}{
Size of the oligonucleotide (e.g.for tetranucleotide it is 4,for hexanucleotide it is 6 ). Default is set to 1.
}
  \item{norm}{
normalization of oligonucleotide frequency (OF) Profile (0 - no normalization, 1 - normalize the OF profile). Default is set to norm = 0.
}
}
\details{
Oligonucleotide Frequency Derived Error Gradient (OFDEG) attempts to capture the convergence behavior by subsampling the genomic fragment and measuring the decrease in error as the length of the subsamples increases upto the fragment lenth.OFDEG, derived from the oligonucleotide frequency profile of a DNA sequence shows that it is possible to obtain a meaningful phylogenetic signal for relatively short DNA sequences.
}
\value{
This function returns a data frame containing error gradients of each nucleotide sequence.
}
\references{
Saeed, I., Halgamuge, S.K. The oligonucleotide frequency derived error gradient and its application to the binning of metagenome fragments. BMC Genomics 10, S10 (2009).
}
\author{
Dr. Anu Sharma,
Dr. Sanjeev Kumar
}
\examples{
library(GenomicSig)
\donttest{
OFDEG(sequence= system.file("extdata/Nuc_sequence.fasta", package = "GenomicSig")[1], c=60, rc=10
, d=10, m=50, t=6, k=1, norm=0)
}
}
