requireNamespace("seqinr", "Biostrings")
RSCU <- function(sequence) {
  # Function to calculate RSCU for a single sequence
  RSCU1 <- function(data) {
    rsu <- function(seq) {
      y <- unlist(strsplit(seq, ""))  # Split the sequence into individual characters
      z <- seqinr::uco(y, index = c("rscu"))  # Assuming uco() is a function to calculate RSCU
      z <- round(z, 4)  # Round to 4 decimal places
      return(z)
    }

    f_res <- lapply(data, rsu)  # Apply RSCU function to each sequence
    d <- data.frame(f_res)  # Convert list to data frame
    e <- t(d)  # Transpose the data frame
    return(e)
  }

  # Read the FASTA file
  fasta_sequences <- Biostrings::readDNAStringSet(sequence)
  sequence_data <- as.character(fasta_sequences)  # Convert to character vector

  # Generate RSCU content
  res <- RSCU1(sequence_data)
  return(res)
}
