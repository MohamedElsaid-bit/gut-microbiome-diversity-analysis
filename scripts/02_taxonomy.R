#!/usr/bin/env Rscript
# Assign genus level taxonomy to the ASVs with DADA2's naive Bayesian
# classifier against the real SILVA v138.1 reference training set.
#
# Run from the repo root, after 01_dada2_asv_inference.R:
#   Rscript scripts/02_taxonomy.R

suppressMessages(library(dada2))

seqtab.nochim <- readRDS("data/processed/seqtab_nochim.rds")
ref_path <- "data/reference/silva_nr99_v138.1_train_set.fa.gz"

if (!file.exists(ref_path)) {
  stop(
    "SILVA reference not found at ", ref_path,
    ". Download it first, see README 'How to run it'."
  )
}

set.seed(100) # assignTaxonomy's bootstrap step is stochastic
taxa <- assignTaxonomy(seqtab.nochim, ref_path, multithread = 4)

saveRDS(taxa, "data/processed/taxonomy.rds")

taxa_print <- taxa
rownames(taxa_print) <- NULL
write.csv(taxa_print, "results/tables/asv_taxonomy.csv")

n_unclassified_phylum <- sum(is.na(taxa[, "Phylum"]))
cat(sprintf(
  "Assigned taxonomy to %d ASVs; %d unclassified at Phylum level\n",
  nrow(taxa), n_unclassified_phylum
))
cat("\nPhylum level composition (number of ASVs):\n")
print(table(taxa[, "Phylum"], useNA = "ifany"))
