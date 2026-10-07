#!/usr/bin/env Rscript
# ASV inference with DADA2 on the real mouse gut (fecal) 16S V4 amplicon
# reads from Kozich et al. 2013 (the "MiSeq SOP" dataset), 9 early and
# 10 late post-weaning samples. The Mock community sample is excluded
# here; it is a sequencing control, not a biological sample.
#
# Run from the repo root: Rscript scripts/01_dada2_asv_inference.R

suppressMessages(library(dada2))

raw_dir <- "data/raw/MiSeq_SOP"
filt_dir <- "data/processed/filtered"
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

fnFs <- sort(list.files(raw_dir, pattern = "_R1_001.fastq$", full.names = TRUE))
fnRs <- sort(list.files(raw_dir, pattern = "_R2_001.fastq$", full.names = TRUE))

# Exclude the Mock community control (not a biological sample).
keep <- !grepl("^Mock", basename(fnFs))
fnFs <- fnFs[keep]
fnRs <- fnRs[keep]

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
stopifnot(length(fnFs) == 19, length(fnRs) == 19)

# Quality profiles, inspected before picking truncation lengths below.
png("results/figures/quality_profile_forward.png", width = 1000, height = 700, res = 120)
print(plotQualityProfile(fnFs[1:4]))
dev.off()
png("results/figures/quality_profile_reverse.png", width = 1000, height = 700, res = 120)
print(plotQualityProfile(fnRs[1:4]))
dev.off()

filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# truncLen chosen from the quality profiles above: forward quality stays
# high through ~240 bp; reverse quality (as is typical for Illumina
# paired-end runs) degrades earlier, so it is truncated harder at 160 bp.
# This still leaves > 20 bp of overlap for the ~253 bp V4 amplicon.
cat("STEP: filterAndTrim\n"); flush(stdout())
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(240, 160),
  maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = 4
)
write.csv(out, "results/tables/filter_trim_summary.csv")

cat("STEP: learnErrors forward\n"); flush(stdout())
errF <- learnErrors(filtFs, multithread = 4)
saveRDS(errF, "data/processed/errF.rds")
cat("STEP: learnErrors reverse\n"); flush(stdout())
errR <- learnErrors(filtRs, multithread = 4)
saveRDS(errR, "data/processed/errR.rds")
png("results/figures/error_rates_forward.png", width = 1000, height = 800, res = 120)
print(plotErrors(errF, nominalQ = TRUE))
dev.off()

cat("STEP: dada forward\n"); flush(stdout())
dadaFs <- dada(filtFs, err = errF, multithread = 4)
cat("STEP: dada reverse\n"); flush(stdout())
dadaRs <- dada(filtRs, err = errR, multithread = 4)

cat("STEP: mergePairs\n"); flush(stdout())
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "data/processed/seqtab_raw.rds")
cat(sprintf("ASV table: %d samples x %d ASVs\n", nrow(seqtab), ncol(seqtab)))
cat("ASV length distribution:\n")
print(table(nchar(getSequences(seqtab))))

cat("STEP: removeBimeraDenovo\n"); flush(stdout())
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 4, verbose = TRUE)
cat(sprintf(
  "After chimera removal: %d ASVs kept (%.1f%% of reads retained)\n",
  ncol(seqtab.nochim), 100 * sum(seqtab.nochim) / sum(seqtab)
))

getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs, getN), sapply(dadaRs, getN),
  sapply(mergers, getN), rowSums(seqtab.nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, "results/tables/reads_tracked_through_pipeline.csv")

saveRDS(seqtab.nochim, "data/processed/seqtab_nochim.rds")
cat("\nWrote data/processed/seqtab_nochim.rds and results/tables/reads_tracked_through_pipeline.csv\n")
