library(ggplot2)
library(dplyr)

evec_file <- "pca_no_Osillanus.eigenvec"
eval_file <- "pca_no_Osillanus.eigenval"
popmap_file <- "config/popmap.tsv"

pca <- read.table(evec_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pca)[1:2] <- c("FID", "IID")
colnames(pca)[3:ncol(pca)] <- paste0("PC", seq_len(ncol(pca) - 2))

pca <- pca %>%
  mutate(PC1 = as.numeric(PC1), PC2 = as.numeric(PC2))

eig <- scan(eval_file, quiet = TRUE)
var_exp <- eig / sum(eig) * 100

popmap <- read.table(popmap_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(popmap) <- c("IID", "POP")

df <- pca %>% left_join(popmap, by = "IID")

pop_levels <- c(
  "GSDW", "IJWH", "IJYD", "IJBD", "IJHG", "SOSA", "MTGA",
  "IJBP", "YYOS", "MTSA", "CCGE", "IJGS", "HCGW", "HCGU",
  "MTCA", "SCDI", "MTGY"
)

df <- df %>% mutate(POP = factor(POP, levels = pop_levels))

pop_cols <- c(
  GSDW = "#08306B", IJWH = "#2171B5", IJYD = "#6BAED6", IJBD = "#BDD7E7",
  IJHG = "#3F007D", SOSA = "#54278F", MTGA = "#6A51A3", IJBP = "#807DBA",
  YYOS = "#9E9AC8", MTSA = "#BCBDDC", CCGE = "#DADAEB", IJGS = "#756BB1",
  HCGW = "#8856A7", HCGU = "#8C6BB1", MTCA = "#B358A3",
  SCDI = "#CB181D", MTGY = "#FB6A4A"
)

hull_df <- df %>%
  group_by(POP) %>%
  filter(n() >= 3) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

p <- ggplot(df, aes(PC1, PC2, color = POP)) +
  geom_polygon(
    data = hull_df,
    aes(fill = POP, group = POP),
    alpha = 0.15,
    color = NA,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(size = 4.8) +
  scale_color_manual(values = pop_cols, breaks = pop_levels, limits = pop_levels) +
  scale_fill_manual(values = pop_cols, breaks = pop_levels, limits = pop_levels) +
  labs(
    x = sprintf("PC1 (%.2f%%)", var_exp[1]),
    y = sprintf("PC2 (%.2f%%)", var_exp[2]),
    color = NULL
  ) +
  theme_bw(base_size = 17, base_family = "Times New Roman") +
  theme(
    axis.title = element_text(size = 21, color = "black"),
    axis.text = element_text(size = 18, color = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 17, color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.height = grid::unit(0.75, "cm"),
    legend.spacing.y = grid::unit(0.12, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5.5)))

ggsave("PCA_no_Osillanus_publication.pdf", p, width = 13, height = 8.2, device = cairo_pdf)
ggsave("PCA_no_Osillanus_publication.png", p, width = 13, height = 8.2, dpi = 300)
