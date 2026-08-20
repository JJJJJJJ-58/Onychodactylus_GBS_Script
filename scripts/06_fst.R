#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(hierfstat))

raw_file <- "gbs_pruned.raw"
pop_file <- "config/popmap.tsv"
nperm <- 999
outpref <- "pairwise_fst"

raw <- read.table(raw_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
iid <- raw$IID
geno <- raw[, -(1:6)]
geno[geno == -9] <- NA

popmap <- read.table(pop_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(popmap) <- c("IID", "POP")
m <- match(iid, popmap$IID)
stopifnot(!any(is.na(m)))
grp <- factor(popmap$POP[m])

g <- as.matrix(geno)
storage.mode(g) <- "numeric"
g[g == 0] <- 11
g[g == 1] <- 12
g[g == 2] <- 22

hf <- data.frame(pop = as.integer(grp), g, check.names = FALSE)
fst_obs <- pairwise.WCfst(hf)

set.seed(1)
K <- nrow(fst_obs)
cnt <- matrix(0, nrow = K, ncol = K, dimnames = dimnames(fst_obs))
is_na <- is.na(fst_obs)
pop0 <- hf$pop

for (b in seq_len(nperm)) {
  hf$pop <- sample(pop0, replace = FALSE)
  fst_b <- pairwise.WCfst(hf)
  ok <- !is_na & !is.na(fst_b)
  cnt[ok] <- cnt[ok] + (fst_b[ok] >= fst_obs[ok])
}

hf$pop <- pop0
p_mat <- (cnt + 1) / (nperm + 1)
diag(p_mat) <- NA

p_vec <- p_mat[upper.tri(p_mat)]
p_adj <- p.adjust(p_vec, method = "BH")
p_adj_mat <- matrix(NA_real_, K, K, dimnames = dimnames(p_mat))
p_adj_mat[upper.tri(p_adj_mat)] <- p_adj
p_adj_mat <- t(p_adj_mat)
p_adj_mat[upper.tri(p_adj_mat)] <- p_adj

write.table(fst_obs, paste0(outpref, ".tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(p_mat, paste0(outpref, "_p_perm.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(p_adj_mat, paste0(outpref, "_p_perm_BH.tsv"), sep = "\t", quote = FALSE, col.names = NA)
