#!/usr/bin/env Rscript
# Alpha diversity, beta diversity (PCoA + PERMANOVA), and a phylum level
# taxonomy bar chart, comparing the real Early (days 0-9 post weaning)
# vs Late (days 141-150) mouse gut samples.
#
# Run from the repo root, after 02_taxonomy.R:
#   Rscript scripts/03_phyloseq_diversity.R

suppressMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
})

seqtab.nochim <- readRDS("data/processed/seqtab_nochim.rds")
taxa <- readRDS("data/processed/taxonomy.rds")
metadata <- read.table("config/metadata.tsv", header = TRUE, sep = "\t", row.names = "sample_id")

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(metadata),
  tax_table(taxa)
)
# Replace the long DNA sequences used as ASV names with short IDs, keeping
# the sequences themselves as a reference table rather than row names.
asv_seqs <- taxa_names(ps)
asv_ids <- paste0("ASV", seq_along(asv_seqs))
write.csv(data.frame(asv_id = asv_ids, sequence = asv_seqs), "results/tables/asv_sequences.csv", row.names = FALSE)
taxa_names(ps) <- asv_ids

cat(sprintf("phyloseq object: %d samples, %d ASVs\n", nsamples(ps), ntaxa(ps)))

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# --- Alpha diversity: Shannon and Chao1, Early vs Late -----------------
alpha <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
alpha$sample_id <- rownames(alpha)
alpha$time_group <- sample_data(ps)$time_group[match(alpha$sample_id, rownames(sample_data(ps)))]
write.csv(alpha, "results/tables/alpha_diversity.csv", row.names = FALSE)

shannon_test <- wilcox.test(Shannon ~ time_group, data = alpha)
chao1_test <- wilcox.test(Chao1 ~ time_group, data = alpha)
alpha_tests <- data.frame(
  metric = c("Shannon", "Chao1"),
  test = "Wilcoxon rank-sum",
  W = c(shannon_test$statistic, chao1_test$statistic),
  p_value = c(shannon_test$p.value, chao1_test$p.value)
)
write.csv(alpha_tests, "results/tables/alpha_diversity_tests.csv", row.names = FALSE)
cat("\nAlpha diversity, Early vs Late (Wilcoxon):\n")
print(alpha_tests)

p_alpha <- ggplot(alpha, aes(x = time_group, y = Shannon, fill = time_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, height = 0, size = 2) +
  scale_fill_manual(values = c("Early" = "#2a78d6", "Late" = "#eb6834")) +
  labs(title = "Shannon alpha diversity: early vs late",
       x = NULL, y = "Shannon index") +
  theme_bw(base_size = 13) + theme(legend.position = "none")
ggsave("results/figures/alpha_diversity_shannon.png", p_alpha, width = 5.5, height = 4.5, dpi = 150)

p_chao1 <- ggplot(alpha, aes(x = time_group, y = Chao1, fill = time_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, height = 0, size = 2) +
  scale_fill_manual(values = c("Early" = "#2a78d6", "Late" = "#eb6834")) +
  labs(title = "Chao1 richness: early vs late",
       x = NULL, y = "Chao1 estimated richness") +
  theme_bw(base_size = 13) + theme(legend.position = "none")
ggsave("results/figures/alpha_diversity_chao1.png", p_chao1, width = 5.5, height = 4.5, dpi = 150)

# --- Beta diversity: Bray-Curtis, PCoA, PERMANOVA ----------------------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
bray_dist <- phyloseq::distance(ps_rel, method = "bray")

set.seed(42)
permanova <- adonis2(bray_dist ~ time_group, data = as(sample_data(ps), "data.frame"), permutations = 999)
write.csv(as.data.frame(permanova), "results/tables/permanova_bray_curtis.csv")
cat("\nPERMANOVA (Bray-Curtis ~ time_group, 999 permutations):\n")
print(permanova)

ord <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)
eig_pct <- 100 * ord$values$Relative_eig[1:2]
ord_df <- data.frame(
  Axis1 = ord$vectors[, 1], Axis2 = ord$vectors[, 2],
  time_group = sample_data(ps)$time_group
)
p_pcoa <- ggplot(ord_df, aes(x = Axis1, y = Axis2, color = time_group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Early" = "#2a78d6", "Late" = "#eb6834")) +
  labs(
    title = "PCoA of Bray-Curtis dissimilarity",
    subtitle = sprintf("PERMANOVA R2 = %.3f, p = %.3f (999 permutations)",
                        permanova$R2[1], permanova$`Pr(>F)`[1]),
    x = sprintf("Axis 1 (%.1f%%)", eig_pct[1]),
    y = sprintf("Axis 2 (%.1f%%)", eig_pct[2]),
    color = "Time group"
  ) +
  theme_bw(base_size = 13)
ggsave("results/figures/pcoa_bray_curtis.png", p_pcoa, width = 6, height = 5, dpi = 150)

# --- Phylum level stacked taxonomy bar chart ---------------------------
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
phylum_df <- psmelt(ps_phylum)
phylum_df$Phylum <- as.character(phylum_df$Phylum)
phylum_df$Phylum[is.na(phylum_df$Phylum)] <- "Unclassified"

# Collapse low-abundance phyla into "Other" for a readable legend.
mean_abund <- aggregate(Abundance ~ Phylum, phylum_df, mean)
top_phyla <- mean_abund$Phylum[order(-mean_abund$Abundance)][1:6]
phylum_df$Phylum <- ifelse(phylum_df$Phylum %in% top_phyla, phylum_df$Phylum, "Other")

p_tax <- ggplot(phylum_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(~time_group, scales = "free_x", space = "free_x") +
  labs(title = "Phylum level relative abundance, early vs late",
       x = NULL, y = "Relative abundance") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave("results/figures/taxonomy_barchart_phylum.png", p_tax, width = 9, height = 5.5, dpi = 150)

saveRDS(ps, "data/processed/phyloseq_object.rds")
cat("\nWrote alpha/beta diversity tables and figures to results/\n")
