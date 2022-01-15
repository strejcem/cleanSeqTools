### Usage:
### bash: conda create -n its cutadapt r-base bioconductor-dada2 bioconductor-shortread bioconductor-biostrings
### bash: conda activate its
### bash: Rscript cutITSprimers.R | tee cut_ITS_primers.log

## Michal Strejcek @ UCT Prague
## modified from https://benjjneb.github.io/dada2/ITS_workflow.html
## v1.1 March 22, 2021

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

path <- "renamed"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

FWD <- "AACTTTYRRCAAYGGATCWCT"
REV <- "AGCCTCCGCTTATTGATATGCTTAART" 

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path("tmp", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path("tmp", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#pick random sample
sample_pick <- sample(length(fnFs.filtN), 1)
fnFs.filtN[[sample_pick]]

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[sample_pick]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[sample_pick]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[sample_pick]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[sample_pick]]))

cutadapt <- "cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path("removed_primers")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste0("-a ^", FWD, "...", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste0("-A ^", REV, "...", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "--discard-untrimmed",
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[100]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[100]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[100]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[100]]))

message("Removing temporary files")
unlink("fnFs.filtN", recursive=TRUE)
unlink("fnRs.filtN", recursive=TRUE)