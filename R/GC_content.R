requireNamespace("Biostrings", "seqinr")
GC_content <- function(sequence) {

  GC_content1 <- function(data) {

    calculate_gc <- function(seq) {
      # Split the sequence into individual characters
      char_seq <- unlist(strsplit(seq, split = ""))
      # Calculate GC content
      res <- seqinr::GC(char_seq)
      val <- round(res, 4)
      return(val)
    }

    gc_results <- lapply(data, calculate_gc)
    gc_df <- data.frame(GC_content = unlist(gc_results))
    # rownames(gc_df) <- NULL  # Remove row names
    return(gc_df)
  }

  # Read the FASTA file
  fasta_sequences <- Biostrings::readDNAStringSet(sequence)
  sequence_data <- as.character(fasta_sequences)

  # Generate GC content
  res <- GC_content1(sequence_data)
  return(res)
}
