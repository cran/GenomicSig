requireNamespace("Biostrings", "stats")

OFDEG <- function(sequence, c, rc, d, m, t, k, norm=0) {

  OFDEG1 <- function(data, c, rc, d, m, t, k, norm) {

    extract_subseq <- function(data, j, c, rc) {
      seq <- as.character(data[[j]])
      L <- nchar(seq)  # Use nchar to get the length of the sequence

      # Ensure sampling is valid
      if (L - c + 1 < 1) {
        return(NULL)  # Return NULL if not enough sequence length
      }

      # Adjust the sample size if necessary
      size_to_sample <- min(L - c + 1, rc)
      ind <- sample(L - c + 1, size_to_sample)  # Adjust sample size

      subseq_list <- lapply(ind, function(x) {
        substr(seq, x, x + c - 1)
      })

      return(subseq_list)
    }

    ofdeg_fun <- function(subseq, c, d, m, t, k, norm = 0) {

      mean_error <- function(x, subseq, c, d, k, norm) {
        L <- nchar(subseq)  # Get the length of the subsequence

        if (L - x + 1 >= d) {
          ind1 <- sample(L - x + 1, d)
        } else if (L - x + 1 > 0) {
          ind1 <- seq_len(L - x + 1)  # Adjusted to sample all if population is smaller
        } else {
          return(NA)  # Return NA if there's no valid subsequence length
        }

        subseq_list1 <- lapply(ind1, function(i) {
          substr(subseq, i, i + x - 1)
        })

        subseq_Sc <- Biostrings::DNAStringSet(paste(subseq, collapse = ""))
        ofp_Sc <- Biostrings::oligonucleotideFrequency(subseq_Sc, k)
        ofp_Sc <- if (norm == 1) ofp_Sc / sum(ofp_Sc) else ofp_Sc

        error_fn <- function(x) {
          seq_Si <- Biostrings::DNAStringSet(paste(x, collapse = ""))
          ofp_Si <- Biostrings::oligonucleotideFrequency(seq_Si, k)
          ofp_Si <- if (norm == 1) ofp_Si / sum(ofp_Si) else ofp_Si
          sqrt(sum((ofp_Sc - ofp_Si)^2))
        }

        err_d <- sapply(subseq_list1, error_fn)
        mean(err_d, na.rm = TRUE)  # Calculate mean, ignoring NAs
      }

      m_ind <- seq.int(from = m, to = c, by = t)
      mean_error_ind <- sapply(m_ind, mean_error, subseq = subseq, c = c, d = d, k = k, norm = norm)
      ofdeg <- stats::lm(mean_error_ind ~ m_ind)

      list(ofdeg$coefficients["m_ind"], m_ind, mean_error_ind)
    }

    ofdeg_rob_est <- vector("numeric", length(data))

    for (j in seq_along(data)) {
      subseq_all <- extract_subseq(data, j, c, rc)

      if (is.null(subseq_all)) {
        next  # Skip if no valid subsequence is available
      }

      ofdeg_all <- numeric(rc)
      m_all <- error_all <- list()

      for (i in seq_len(rc)) {
        z <- ofdeg_fun(subseq_all[[i]], c, d, m, t, k, norm)
        ofdeg_all[i] <- z[[1]]
        m_all[[i]] <- z[[2]]
        error_all[[i]] <- z[[3]]
      }

      error_all <- unlist(error_all)
      m_all <- unlist(m_all)
      ofdeg_robust <- stats::lm(error_all ~ m_all)
      ofdeg_rob_est[j] <- ofdeg_robust$coefficients["m_all"]
    }

    OFDEG_Value <- data.frame(Ofdeg_Robust = unlist(ofdeg_rob_est))
    rownames(OFDEG_Value) <- toupper(names(data))
    colnames(OFDEG_Value) <- toupper(colnames(OFDEG_Value))

    return(OFDEG_Value)
  }

  # Read the FASTA file
  fasta_sequences <- Biostrings::readDNAStringSet(sequence)
  sequence_data <- as.character(fasta_sequences)

  # Call OFDEG1 with appropriate parameters
  res <- OFDEG1(sequence_data, c, rc, d, m, t, k, norm)
  return(res)
}
