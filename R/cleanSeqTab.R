message("cleanSeqTools.R v1.28 loaded (16-07-2020)")
message("Depencencies: 'dada2', 'ShortRead', 'phyloseq'")

###############################################################################
#******************************************************************************
###############################################################################

renameByID <- function(dir.from,
                       dir.to,
                       extension = "fastq.gz",
                       csvFile = NULL,
                       ID.from = NULL,
                       ID.to = NULL,
                       prefix = NULL) {
  if (is.null(c(csvFile,
                ID.from,
                ID.to))) {
    stop("Either 'csvFile' or 'ID.from' and 'ID.to' must be provided!")
  }
  
  if (is.null(csvFile)) {
    if (is.null(ID.from) | is.null(ID.to)) {
      stop("Both 'ID.from' and 'ID.to' arguments must be provided!")
    }
  }
  
  if (!is.null(csvFile)) {
    sampleTab <- read.csv(csvFile,
                          row.names = 1,
                          stringsAsFactors = F)
    
    if (!all(dim(sampleTab)) > 0) {
      stop("Something wrong with the '", basename(csvFile), "'")
    }
    
    samples_new <- unname(unlist(sampleTab))
    samples_new[samples_new == ""] <- NA
    notEmpty <- !is.na(samples_new)
    message("Expecting ",
            sum(notEmpty),
            " samples described in ",
            basename(csvFile))
    x <- dim(sampleTab)[1]
    y <- dim(sampleTab)[2]
    samples_old <-
      paste0(rep(LETTERS[1:x], times = y), rep(sprintf("%02d", 1:y), each = x))
    samples_old <- samples_old[notEmpty]
    samples_new <- samples_new[notEmpty]
  }
  
  if (!(is.null(ID.from) & is.null(ID.to))) {
    samples_old <- ID.from
    samples_new <- ID.to
  }
  
  if (any(duplicated(samples_new))) {
    stop("New sample names must be unique!\n", samples_new[duplicated(samples_new)])
  }
  
  fln <-
    list.files(path = dir.from,
               pattern = extension,
               full.names = T)
  
  if (!grepl("_R[1|2]", fln[2])) {
    stop("Sequence files are expected to have _R1 and _R2 identifier")
  }
  
  for (i in seq_along(samples_old)) {
    files <- fln[grepl(samples_old[i], fln)]
    if (length(files) != 2) {
      stop(
        "Irregular number of samples detected, should be exactly 2 (R1 and R2):\n",
        paste(files, collapse = ", ")
      )
    }
    suffix <- sub("^.*(_R[1|2]\\..*$)", "\\1", files)
    file.rename(from = file.path(dir.from, files),
                to = file.path(dir.to, paste0(prefix, samples_new[i], suffix)))
  }
  message("Done")
}

###############################################################################
#******************************************************************************
###############################################################################

checkPrimers <-
  function(dir,
           fPrimer = "GTGYCAGCMGCNGCGG",
           rPrimer = "CCGYCAATTYMTTTRAGTTT") {
    fl <- list.files(dir)
    
    rndSample <- sample(length(fl) / 2, 1)
    
    for (i in 1:2) {
      fls <-
        list.files(dir,
                   pattern = paste0("_R", i, "\\."),
                   full.names = T)
      
      message(paste0("Randomly choosen: ", basename(fls[rndSample])))
      
      fq <- ShortRead::readFastq(dirPath = fls[rndSample])
      reads <- ShortRead::sread(fq)
      
      primer <- ifelse(i == 1, fPrimer, rPrimer)
      prim <- ifelse(i == 1, "forward", "reverse")
      
      message(paste0(
        "Checking ",
        prim,
        " primer. Expected peak at ",
        nchar(primer),
        "bp"
      ))
      
      primerMatch <- Biostrings::trimLRPatterns(primer,
                                                subject = reads,
                                                Lfixed = F,
                                                ranges = T)
      print(table(primerMatch@start - 1))
      
      message("\nTop 10 20-mers:")
      first20 <- substr(reads, start = 1, stop = 20)
      tbl <- sort(table(first20), decreasing = T)
      df <- data.frame(tbl)
      print(df[1:10,])
    }
  }

###############################################################################
#******************************************************************************
###############################################################################

removePrimers <-
  function(dir.from,
           dir.to,
           fPrimer = "GTGYCAGCMGCNGCGG",
           rPrimer = "CCGYCAATTYMTTTRAGTTT",
           numDiff = 0,
           parallel = FALSE) {
    fls <- list.files(path = dir.from,
                      ignore.case = T,
                      full.names = T)
    if (!dir.exists(dir.to)) {
      dir.create(dir.to)
    }
    
    sample.rename <-
      function(x, dir.to, fPrimer, rPrimer, numDiff) {
        output <- NULL
        sampleName <- basename(x)
        f <- ShortRead::readFastq(dirPath = x)
        reads.in <- length(f)
        if (reads.in == 0) {
          message("Empty sample (discarded): ", sampleName)
          break
        }
        
        toTrim <- if (grepl("_R1\\.", sampleName)) {
          as.character(fPrimer)
        } else if (grepl("_R2\\.", x)) {
          as.character(rPrimer)
        } else
          stop(paste0("No 'R1' or 'R2' identifier in ", x))
        
        range <- Biostrings::trimLRPatterns(
          Lpattern = toTrim,
          subject = f,
          Lfixed = FALSE,
          ranges = TRUE,
          max.Lmismatch = rep(numDiff, nchar(toTrim))
        )
        
        f.new <- f[range@start == nchar(toTrim) + 1]
        f.new <- Biostrings::trimLRPatterns(
          Lpattern = toTrim,
          subject = f.new,
          Lfixed = FALSE,
          max.Lmismatch = rep(numDiff, nchar(toTrim))
        )
        removed <- sum(range@start != nchar(toTrim) + 1)
        percent <- round(removed * 100 / reads.in, 1)
        
        ## TODO: check if all(writeFasta == files)
        invisible(ShortRead::writeFastq(
          f.new,
          file = file.path(dir.to, sampleName),
          compress = grepl("gz$", x)
        ))
        output <-
          rbind(
            output,
            data.frame(
              "sample" = sampleName,
              "barPrimerLen" = nchar(toTrim),
              "raw.reads" = reads.in,
              "removed" = removed,
              "removedPerc" = percent
            )
          )
        return(output)
      }
    if (is.numeric(parallel)) {
      cl <- parallel::makeCluster(parallel)
      parallel <- TRUE
    } else if (isTRUE(parallel)) {
      cl <- parallel::makeCluster(parallel::detectCores() - 1)
    } else if (!isFALSE(parallel)) {
      stop("Required 'parallel argument values: number of cores, TRUE, FALSE")
    }
    
    if (parallel) {
      output1 <-
        parallel::parLapply(cl, fls, sample.rename, dir.to, fPrimer, rPrimer, numDiff)
      parallel::stopCluster(cl)
    } else {
      output1 <-
        lapply(fls, sample.rename, dir.to, fPrimer, rPrimer, numDiff)
    }
    return(do.call(rbind, output1))
  }

###############################################################################
#******************************************************************************
###############################################################################

cleanSeqTab <-
  function (seqTab,
            numDiffs = 1,
            parentRatio = 1,
            parallel = TRUE,
            KmerNWdist = NULL,
            kmer_size = 5) {
    if (!is.matrix(seqTab))
      stop(paste0("Object '", deparse(substitute(seqTab)),
                  "' is not a matrix!"))
    if (numDiffs > 2)
      message(
        "'numDiffs' was set to be higher than 2!\n This function is not optimized for such high value!"
      )
    if (is.null(KmerNWdist)) {
      seqs <- colnames(seqTab)
      message("Calculating Kmer-NW distance")
      KmerNWdist <- hammingFast(seqs, kmer_size, numDiffs,
                                parallel)
    }
    else
      message("Applying supplemented Kmer-NW distance matrix!")
    KmerNWdist <- KmerNWdist[KmerNWdist$diff_bases <= numDiffs,]
    abundance <- colSums(seqTab)
    sequenceIdxToProcess <- unique(c(KmerNWdist$s1, KmerNWdist$s2))
    message("Sequence table cleaning")
    
    while (length(sequenceIdxToProcess) > 0) {
      lowestAbundanceIdx <- which.min(abundance[sequenceIdxToProcess])
      IdxToProcess <- sequenceIdxToProcess[lowestAbundanceIdx]
      seqPairs <- c(KmerNWdist[KmerNWdist$s2 == IdxToProcess,
                               "s1"], KmerNWdist[KmerNWdist$s1 == IdxToProcess,
                                                 "s2"])
      if (length(seqPairs) == 0) {
        sequenceIdxToProcess <- sequenceIdxToProcess[-lowestAbundanceIdx]
        next
      }
      abundance <- colSums(seqTab)
      parent_max <- seqPairs[which.max(abundance[seqPairs])]
      seqPairsMessage <- paste(seqPairs, abundance[seqPairs],
                               sep = "(")
      seqPairsMessage <-
        paste0(seqPairsMessage, ")", collapse = ", ")
      if (abundance[parent_max] / abundance[IdxToProcess] >=
          parentRatio) {
        seqTab[, parent_max] <- seqTab[, parent_max] + seqTab[,
                                                              IdxToProcess]
        seqTab[, IdxToProcess] <- 0
        message(
          paste0(
            "Sequence #",
            IdxToProcess,
            "(",
            abundance[IdxToProcess],
            ")",
            " was matched to:",
            seqPairsMessage,
            " and was merged with #",
            parent_max
          )
        )
      }
      sequenceIdxToProcess <-
        sequenceIdxToProcess[-lowestAbundanceIdx]
    }
    message("Done!")
    seqTab <- seqTab[, colSums(seqTab) > 0]
    seqTab <- seqTab[, order(colSums(seqTab), decreasing = T)]
    return(seqTab)
  }




hammingFast <-
  function (seqs,
            kmer_size = 5,
            numDiffs = 1,
            parallel = TRUE) {
    hamming_fun <- function(i, n, seqs, kmer_size, numDiffs) {
      Kmers <- min(nchar(seqs)) - kmer_size + 1
      Kmer_cutoff <- kmer_size * numDiffs / Kmers
      Kmer_cutoff <- min(Kmer_cutoff, 0.5)
      s1 <- rep.int(i, times = n - i)
      s2 <- seq.int(i + 1, n)
      kmer_dist <- dada2:::kmer_dist(seqs[s1], seqs[s2], kmer_size)
      kmer_test <- kmer_dist <= Kmer_cutoff
      if (sum(kmer_test) == 0) {
        return(NULL)
      }
      NW <- dada2:::C_nwvec(
        seqs[s1[kmer_test]],
        seqs[s2[kmer_test]],
        match = 5,
        mismatch = -4,
        gap_p = -8,
        band = 16,
        endsfree = FALSE
      )
      lenNW <- length(NW) / 2
      diff_bases <- sapply(seq_len(lenNW), function(x) {
        xx = unlist(strsplit(NW[2 * x - 1], ""))
        yy = unlist(strsplit(NW[2 * x], ""))
        sum(xx != yy)
      })
      output <- data.frame(
        s1 = s1[kmer_test],
        s2 = s2[kmer_test],
        Kmer_dist = kmer_dist[kmer_test],
        diff_bases = diff_bases
      )
      return(output)
    }
    
    stopifnot(class(seqs) == "character")
    n <- length(seqs)
    
    if (is.numeric(parallel)) {
      cl <- parallel::makeCluster(parallel)
      parallel <- TRUE
    } else if (isTRUE(parallel)) {
      cl <- parallel::makeCluster(parallel::detectCores() - 1)
    } else if (!isFALSE(parallel)) {
      stop("Required 'parallel argument values: number of cores, TRUE, FALSE")
    }
    
    if (parallel) {
      hamming <- parallel::parLapply(cl, seq_len(n - 1), hamming_fun,
                                     n, seqs, kmer_size, numDiffs)
      parallel::stopCluster(cl)
    }
    else {
      hamming <- lapply(seq_len(n - 1), hamming_fun, n, seqs,
                        kmer_size, numDiffs)
    }
    
    hamming <- do.call(rbind, hamming)
    return(hamming)
  }

###############################################################################
#******************************************************************************
###############################################################################

cleanReplicates <-
  function (seqTab,
            replicates,
            prevalence,
            merge = TRUE,
            eraseEmpty = TRUE) {
    if (!is.matrix(as.matrix(seqTab)))
      stop(paste0(
        "Object '",
        deparse(substitute(seqTab)),
        "' is not a matrix or data.frame!"
      ))
    
    replicates <- as.factor(replicates)
    
    if (class(prevalence) != "numeric")
      stop("'prevalence' must be numeric!")
    
    if (length(prevalence) == 1) {
      prevalence <- rep(prevalence, length(levels(replicates)))
    } else if (length(prevalence != length(levels(replicates)))) {
      stop("'prevalence' must be length of 1 or number of final samples")
    }
    
    selection <- seqTab
    selection[selection > 0] <- 1
    if (merge)
      message("Replicates samples merged into: ",
              paste(levels(replicates),
                    collapse = ", "))
    num_replicates <- table(replicates)
    prevalence[prevalence > num_replicates] <-
      num_replicates[prevalence >
                       num_replicates]
    selection <-
      t(sapply(levels(replicates), function(x)
        colSums(selection[replicates ==
                            x, , drop = FALSE])))
    selection <- apply(selection, 2, ">=", prevalence)
    if (merge) {
      output <-
        t(sapply(levels(replicates), function(x)
          colSums(seqTab[replicates ==
                           x, , drop = FALSE])))
      output <- output * selection
    }
    else if (!merge) {
      selection <- selection[replicates,]
      output <- seqTab * selection
    }
    if (eraseEmpty)
      output <- output[, colSums(output) > 0]
    return(output)
  }

###############################################################################
#******************************************************************************
###############################################################################

checkMock <- function(mock.sample, mock.ref, topTaxa = 50) {
  local_fun <- function(seq, mock.ref) {
    NW <- Biostrings::pairwiseAlignment(mock.ref, seq, type = "global")
    nm <- Biostrings::nmismatch(NW)
    data.frame(closest = names(mock.ref[which.min(nm)]) ,
               mismatches = min(nm))
  }
  mock.sample <-
    sort(mock.sample[mock.sample > 0], decreasing = TRUE)
  topTaxa <- min(topTaxa, length(mock.sample))
  mock.sample.subset <- mock.sample[1:topTaxa]
  seqs <- names(mock.sample.subset)
  l <- lapply(seqs, local_fun, mock.ref)
  output <- do.call(rbind, l)
  output[, "abundance"] <- mock.sample.subset
  output[, "fraction"] <-
    round(mock.sample.subset / sum(mock.sample), 4)
  message("Number exact matches ",
          sum(output[, "mismatches"] == 0),
          " out of ",
          length(mock.ref))
  return(output)
}

###############################################################################
#******************************************************************************
###############################################################################

writeDTB <- function(ps, path) {
  stopifnot(!is.null(path))
  
  if (is.null(phyloseq::access(ps, "refseq"))) {
    seq_name <- NULL
    sequence <- phyloseq::taxa_names(ps)
  } else {
    seq_name <- phyloseq::taxa_names(ps)
    sequence <- as.character(phyloseq::refseq(ps))
  }
  
  if (phyloseq::taxa_are_rows(ps)) {
    otu_tab <- phyloseq::otu_table(ps)
  } else {
    otu_tab <- t(phyloseq::otu_table(ps))
  }
  
  tax_tab <- phyloseq::tax_table(ps)
  
  dtb <- data.frame(seq_name, otu_tab, sequence, tax_tab)
  
  write.table(dtb,
              file = path,
              sep = "\t",
              row.names = FALSE)
}
