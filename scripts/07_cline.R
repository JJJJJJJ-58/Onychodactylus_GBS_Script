library(readr)
library(dplyr)
library(sf)
library(ggplot2)
library(MASS)
library(minpack.lm)

plot_font <- "Times New Roman"
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(TimesNewRoman = grDevices::windowsFont("Times New Roman"))
  plot_font <- "TimesNewRoman"
}

qdat <- read_tsv("structure_q_mt.tsv", show_col_types = FALSE)
coord <- read_tsv("coordinates.tsv", show_col_types = FALSE)

coord_sf <- coord %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(32652)

xy <- st_coordinates(coord_sf)
coord_xy <- coord %>% mutate(x = xy[, 1], y = xy[, 2])

start <- coord_xy %>% filter(pop == "GSDW")
end <- coord_xy %>% filter(pop == "SCDI")

vx <- end$x - start$x
vy <- end$y - start$y
vlen <- sqrt(vx^2 + vy^2)

coord_dist <- coord_xy %>%
  mutate(dist_km = ((x - start$x) * vx + (y - start$y) * vy) / vlen / 1000)

cline_dat <- qdat %>%
  left_join(coord_dist %>% select(pop, dist_km), by = "pop") %>%
  mutate(
    qOK = as.numeric(qOK),
    mtOK = case_when(mt == "OK" ~ 1, mt == "NE" ~ 0, TRUE ~ NA_real_)
  ) %>%
  filter(!pop %in% c("MTUM", "YSSS"), mt %in% c("OK", "NE")) %>%
  arrange(dist_km)

eps_q <- 1e-4
cline_dat <- cline_dat %>%
  mutate(
    qOK_clip = pmin(pmax(qOK, eps_q), 1 - eps_q),
    logit_qOK = qlogis(qOK_clip)
  )

m_cline <- nlsLM(
  qOK ~ 1 / (1 + exp(-4 * (dist_km - center) / width)),
  data = cline_dat,
  start = list(center = 55.74, width = 85.21),
  lower = c(center = min(cline_dat$dist_km), width = 0.1),
  upper = c(center = max(cline_dat$dist_km), width = 500),
  control = nls.lm.control(maxiter = 1000)
)

qOK_center <- coef(m_cline)["center"]
qOK_width <- coef(m_cline)["width"]

m_ibd_linear <- lm(logit_qOK ~ dist_km, data = cline_dat)
m_ibd_quadratic <- lm(logit_qOK ~ dist_km + I(dist_km^2), data = cline_dat)

AICc_fun <- function(model) {
  n <- nobs(model)
  k <- attr(logLik(model), "df")
  AIC(model) + (2 * k * (k + 1)) / (n - k - 1)
}

aicc_tbl <- data.frame(
  model = c("Cline_sigmoid", "IBD_linear", "IBD_quadratic"),
  AIC = c(AIC(m_cline), AIC(m_ibd_linear), AIC(m_ibd_quadratic)),
  AICc = c(AICc_fun(m_cline), AICc_fun(m_ibd_linear), AICc_fun(m_ibd_quadratic))
) %>%
  mutate(
    deltaAICc = AICc - min(AICc),
    weight = exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc))
  ) %>%
  arrange(AICc)

write_csv(aicc_tbl, "NEOK_qOK_AICc_final.csv")

mt_site <- cline_dat %>%
  group_by(pop, dist_km) %>%
  summarise(
    mt_OK = sum(mtOK == 1, na.rm = TRUE),
    mt_NE = sum(mtOK == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dist_km)

m_mt <- glm(cbind(mt_OK, mt_NE) ~ dist_km, data = mt_site, family = binomial)

b0 <- coef(m_mt)[1]
b1 <- coef(m_mt)[2]
mt_center <- -b0 / b1
mt_width <- 4 / b1

pred_x <- data.frame(
  dist_km = seq(min(cline_dat$dist_km), max(cline_dat$dist_km), length.out = 500)
)

set.seed(123)
sim_par <- MASS::mvrnorm(5000, mu = coef(m_cline), Sigma = vcov(m_cline)) %>%
  as.data.frame()
names(sim_par) <- names(coef(m_cline))
sim_par <- sim_par %>% filter(width > 0)

sim_pred_mat <- sapply(seq_len(nrow(sim_par)), function(i) {
  1 / (1 + exp(-4 * (pred_x$dist_km - sim_par$center[i]) / sim_par$width[i]))
})

pred_qOK <- pred_x %>%
  mutate(
    Cline_sigmoid = pmin(pmax(predict(m_cline, newdata = pred_x), 0), 1),
    qOK_lwr = pmin(pmax(apply(sim_pred_mat, 1, quantile, probs = 0.025, na.rm = TRUE), 0), 1),
    qOK_upr = pmin(pmax(apply(sim_pred_mat, 1, quantile, probs = 0.975, na.rm = TRUE), 0), 1),
    IBD_linear = pmin(pmax(plogis(predict(m_ibd_linear, newdata = pred_x)), 0), 1),
    IBD_quadratic = pmin(pmax(plogis(predict(m_ibd_quadratic, newdata = pred_x)), 0), 1)
  )

mt_pred <- predict(m_mt, newdata = pred_x, type = "link", se.fit = TRUE)
pred_mt <- pred_x %>%
  mutate(
    mt_cline = plogis(mt_pred$fit),
    mt_lwr = plogis(mt_pred$fit - 1.96 * mt_pred$se.fit),
    mt_upr = plogis(mt_pred$fit + 1.96 * mt_pred$se.fit)
  )

width_bracket_df <- data.frame(
  marker = c("SNP width", "mtDNA width"),
  xmin = c(qOK_center - abs(qOK_width) / 2, mt_center - abs(mt_width) / 2),
  xmax = c(qOK_center + abs(qOK_width) / 2, mt_center + abs(mt_width) / 2),
  y = c(-0.055, -0.095)
) %>%
  mutate(tick_ymin = y - 0.012, tick_ymax = y + 0.012)

p <- ggplot() +
  geom_ribbon(
    data = pred_qOK,
    aes(dist_km, ymin = qOK_lwr, ymax = qOK_upr),
    fill = "#D73027", alpha = 0.13
  ) +
  geom_ribbon(
    data = pred_mt,
    aes(dist_km, ymin = mt_lwr, ymax = mt_upr),
    fill = "#4575B4", alpha = 0.13
  ) +
  geom_point(
    data = cline_dat,
    aes(dist_km, qOK, fill = "SNP"),
    shape = 21, color = "black", stroke = 0.25, size = 3.2, alpha = 0.95,
    position = position_jitter(width = 0.7, height = 0, seed = 123)
  ) +
  geom_point(
    data = cline_dat,
    aes(dist_km, mtOK, fill = "mtDNA Cytb haplotype"),
    shape = 21, color = "black", stroke = 0.25, size = 3.2, alpha = 0.95,
    position = position_jitter(width = 0.7, height = 0.015, seed = 456)
  ) +
  geom_line(data = pred_qOK, aes(dist_km, Cline_sigmoid, color = "SNP cline"), linewidth = 1.45) +
  geom_line(data = pred_qOK, aes(dist_km, IBD_linear, linetype = "IBD linear"), color = "black", linewidth = 0.85) +
  geom_line(data = pred_qOK, aes(dist_km, IBD_quadratic, linetype = "IBD quadratic"), color = "black", linewidth = 0.85) +
  geom_line(data = pred_mt, aes(dist_km, mt_cline, color = "mtDNA cline"), linewidth = 1.45) +
  geom_vline(xintercept = qOK_center, color = "#D73027", linetype = "dashed", linewidth = 0.85) +
  geom_vline(xintercept = mt_center, color = "#4575B4", linetype = "dashed", linewidth = 0.85) +
  geom_segment(
    data = width_bracket_df %>% filter(marker == "SNP width"),
    aes(x = xmin, xend = xmax, y = y, yend = y),
    color = "#D73027", linewidth = 0.9
  ) +
  geom_segment(
    data = width_bracket_df %>% filter(marker == "SNP width"),
    aes(x = xmin, xend = xmin, y = tick_ymin, yend = tick_ymax),
    color = "#D73027", linewidth = 0.9
  ) +
  geom_segment(
    data = width_bracket_df %>% filter(marker == "SNP width"),
    aes(x = xmax, xend = xmax, y = tick_ymin, yend = tick_ymax),
    color = "#D73027", linewidth = 0.9
  ) +
  geom_segment(
    data = width_bracket_df %>% filter(marker == "mtDNA width"),
    aes(x = xmin, xend = xmax, y = y, yend = y),
    color = "#4575B4", linewidth = 0.9
  ) +
  geom_segment(
    data = width_bracket_df %>% filter(marker == "mtDNA width"),
    aes(x = xmin, xend = xmin, y = tick_ymin, yend = tick_ymax),
    color = "#4575B4", linewidth = 0.9
  ) +
  geom_segment(
    data = width_bracket_df %>% filter(marker == "mtDNA width"),
    aes(x = xmax, xend = xmax, y = tick_ymin, yend = tick_ymax),
    color = "#4575B4", linewidth = 0.9
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.02))) +
  scale_fill_manual(values = c("SNP" = "#D73027", "mtDNA Cytb haplotype" = "#4575B4")) +
  scale_color_manual(values = c("SNP cline" = "#D73027", "mtDNA cline" = "#4575B4")) +
  scale_linetype_manual(values = c("IBD linear" = "solid", "IBD quadratic" = "dotdash")) +
  labs(
    x = "Distance projected along the GSDW–SCDI transect (km)",
    y = expression(Q~score~"("*italic("O. koreanus")*"),"~italic("Cytb")~haplotype~frequency),
    fill = NULL, color = NULL, linetype = NULL
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(shape = 21, color = "black", size = 3.5, stroke = 0.25)),
    color = guide_legend(order = 2, override.aes = list(linewidth = 1.4)),
    linetype = guide_legend(order = 3, override.aes = list(color = "black", linewidth = 0.9))
  ) +
  coord_cartesian(ylim = c(-0.12, 1.05), clip = "off") +
  theme_bw(base_size = 14, base_family = plot_font) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 18, l = 10)
  )

ggsave(
  "NEOK_IBD_SNP_mt.pdf", p,
  width = 9, height = 5.8, units = "in", device = grDevices::cairo_pdf
)

ggsave(
  "NEOK_IBD_SNP_mt.png", p,
  width = 9, height = 5.8, units = "in", dpi = 600,
  device = ragg::agg_png, bg = "white"
)
