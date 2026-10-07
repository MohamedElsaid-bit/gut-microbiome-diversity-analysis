#!/usr/bin/env Rscript
# Genus level differential abundance, Early vs Late, with ANCOM-BC2.
# ANCOM-BC corrects for compositionality (relative abundance data cannot
# be compared with a plain t-test or DESeq2 without that correction)
# and returns a bias-corrected log fold change with its own p-values,
# rather than requiring a separate normalization step.
#
# Run from the repo root, after 03_phyloseq_diversity.R:
#   Rscript scripts/04_differential_abundance.R

suppressMessages({
  library(phyloseq)
  library(ANCOMBC)
})

ps <- readRDS("data/processed/phyloseq_object.rds")
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
cat(sprintf("Testing %d genera, Early (n=9) vs Late (n=10)\n", ntaxa(ps_genus)))

set.seed(42)
out <- ancombc2(
  data = ps_genus,
  fix_formula = "time_group",
  p_adj_method = "BH",
  group = "time_group",
  struc_zero = TRUE,
  n_cl = 1
)

res <- out$res
write.csv(res, "results/tables/ancombc_differential_abundance_genus.csv", row.names = FALSE)

lfc_col <- grep("^lfc_time_group", colnames(res), value = TRUE)[1]
q_col <- grep("^q_time_group", colnames(res), value = TRUE)[1]
diff_col <- grep("^diff_time_group", colnames(res), value = TRUE)[1]

sig <- res[res[[diff_col]] == TRUE, c("taxon", lfc_col, q_col)]
sig <- sig[order(sig[[q_col]]), ]

# Attach the actual genus (and phylum) name for each significant ASV-level
# representative, since ANCOM-BC's own output only carries the ASV ID.
tax_lookup <- as.data.frame(tax_table(ps_genus))
tax_lookup$taxon <- rownames(tax_lookup)
sig <- merge(sig, tax_lookup[, c("taxon", "Phylum", "Genus")], by = "taxon")
sig <- sig[order(sig[[q_col]]), ]
write.csv(sig, "results/tables/ancombc_significant_genera.csv", row.names = FALSE)

cat(sprintf(
  "\n%d of %d genera significantly differentially abundant (q < 0.05, BH corrected)\n",
  nrow(sig), nrow(res)
))
print(sig)
