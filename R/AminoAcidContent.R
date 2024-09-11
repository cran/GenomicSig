requireNamespace("Biostrings")

AminoAcidContent <- function(fasta_file, type = c("DNA", "protein")) {
  type = match.arg(type)

  if (type == "DNA"){
    dna <- Biostrings::readDNAStringSet(fasta_file)
    prot <- Biostrings::translate(dna, genetic.code=Biostrings::GENETIC_CODE, no.init.codon=TRUE, if.fuzzy.codon="solve")
    nm <- names(prot)
    prot1 <- as.character(prot)
  }

  if (type == "protein"){
    prot <- Biostrings::readAAStringSet(fasta_file)
    nm <- names(prot)
    prot1 <- as.character(prot)
  }
  protcheck <- function(prot1) {
    AADict <- c(
      "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X", "*"
    )

    all(strsplit(prot1, split = "")[[1]] %in% AADict)
  }

  if (protcheck(prot1) == FALSE) {
    stop("protein has unrecognized amino acid type")
  }

  # 20 Amino Acid Abbrevation Dictionary from
  # https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties

  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X", "*"
  )

  AAC <- function(prot1){summary(
    factor(strsplit(prot1, split = "")[[1]], levels = AADict),
    maxsum = 22
  ) / nchar(prot1)
  }
  res_aac<-lapply(prot1,AAC)
  names(res_aac) <- nm
  res_aac <- as.data.frame(res_aac)
  res_aac <- t(res_aac)
  return(res_aac)
}

