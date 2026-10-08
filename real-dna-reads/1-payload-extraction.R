# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("ShortRead", "Biostrings"))

library(ShortRead)
library(Biostrings)
library(stringr)

wd <- getwd()
fastq_file_path <- file.path(wd, "data-nbt17", "id20.fastq.gz")
fastq_data <- ShortRead::readFastq(fastq_file_path)
reads <- as.character(ShortRead::sread(fastq_data))
rm(fastq_data)

front_primer <- "AGTGCAACAAGTCAATCCGT"
reverse_primer <- "AATTGAATGCTTGCTTGCCG"
front_primer_revcomp <- as.character(reverseComplement(DNAString(front_primer)))
reverse_primer_revcomp <- as.character(reverseComplement(DNAString(reverse_primer)))

# forward <- str_detect(string = reads, pattern = front_primer) & str_detect(string = reads, pattern = reverse_primer)
# revcomp <- str_detect(string = reads, pattern = reverse_primer_revcomp) & str_detect(string = reads, pattern = front_primer_revcomp)

idx_forward_front <- str_locate(reads, pattern = front_primer)
idx_forward_reverse <- str_locate(reads, pattern = reverse_primer)

idx_revcomp_front <- str_locate(reads, pattern = front_primer_revcomp)
idx_revcomp_reverse <- str_locate(reads, pattern = reverse_primer_revcomp)

valid_forward <- !is.na(idx_forward_front[, "start"]) &
  !is.na(idx_forward_reverse[, "start"]) &
  idx_forward_front[, "end"] < idx_forward_reverse[, "start"]

valid_revcomp <- !is.na(idx_revcomp_reverse[, "start"]) &
  !is.na(idx_revcomp_front[, "start"]) &
  idx_revcomp_reverse[, "end"] < idx_revcomp_front[, "start"]

payload_forward <- rep(NA_character_, length(reads))
payload_forward[valid_forward] <- str_sub(
  reads[valid_forward],
  start = idx_forward_front[valid_forward, "end"] + 1,
  end   = idx_forward_reverse[valid_forward, "start"] - 1
)

payload_revcomp <- rep(NA_character_, length(reads))
payload_revcomp[valid_revcomp] <- str_sub(
  reads[valid_revcomp],
  start = idx_revcomp_reverse[valid_revcomp, "end"] + 1,
  end   = idx_revcomp_front[valid_revcomp, "start"] - 1
)

payload_revcomp[valid_revcomp] <- payload_revcomp[valid_revcomp] |> 
  DNAStringSet() |> 
  reverseComplement() |> 
  as.character()

min_length <- 90

ambiguous <- !is.na(payload_forward) & !is.na(payload_revcomp)
print(sum(ambiguous))
payload_forward[ambiguous] <- NA
payload_revcomp[ambiguous] <- NA

payload <- dplyr::coalesce(payload_forward, payload_revcomp)

payload <- payload[!is.na(payload) & (nchar(payload) >= min_length)]

writeLines(text = payload, con = "data-nbt17/nbt17-reads-preprocessed.txt")

refs <- readLines(con = gzfile("data-nbt17/id20.refs.txt.gz")) |> 
  str_sub(start = 24, end = 133)
table(nchar(refs))
writeLines(text = refs, con = "data-nbt17/nbt17-refs-preprocessed.txt")


