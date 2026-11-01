
# ============================================================
# Example: Combined Manhattan and Miami Plot
# ============================================================
#
# This example generates a combined figure consisting of:
#
#   1. Manhattan plot of csmGmm replicated SNPs
#   2. Miami plot comparing African and European GWAS results
#
# Dataset 1: African females (GCST90551892)
# Dataset 2: European females (iCOGS, GCST010098)
#
# Both datasets use genome build 37.
#
# The replication analysis must be completed before running
# this script.
#
# ============================================================


# 1. Load required packages
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)


# 2. Read GWAS summary statistics
# ------------------------------------------------------------

iCOGS <- read_tsv(
  "data/GCST010098.tsv",
  col_select = c(
    variant_id,
    chromosome,
    base_pair_location,
    beta,
    standard_error,
    effect_allele,
    other_allele,
    p_value
  ),
  progress = FALSE
)

Africa <- read_tsv(
  "data/GCST90551892.tsv.gz",
  col_select = c(
    variant_id,
    chromosome,
    base_pair_location,
    beta,
    standard_error,
    effect_allele,
    other_allele,
    p_value
  ),
  progress = FALSE
)


# 3. Read significant csmGmm SNP results
# ------------------------------------------------------------

sig_snps <- read_csv(
  "results/csmgmm_rep_results.csv",
  show_col_types = FALSE
)

plotData <- sig_snps %>%
  transmute(
    Chr = as.numeric(Chr),
    BP = as.numeric(BP),
    cum_avg_lfdr = as.numeric(cum_avg_lfdr)
  ) %>%
  filter(
    !is.na(Chr),
    !is.na(BP),
    !is.na(cum_avg_lfdr),
    Chr >= 1,
    Chr <= 22,
    BP >= 0,
    cum_avg_lfdr > 0
  )


# 4. Prepare African GWAS data
# ------------------------------------------------------------

Africa_miami <- Africa %>%
  transmute(
    Chr = as.numeric(chromosome),
    BP = as.numeric(base_pair_location),
    logp = -log10(p_value)
  ) %>%
  filter(
    !is.na(Chr),
    !is.na(BP),
    Chr >= 1,
    Chr <= 22,
    BP >= 0,
    is.finite(logp)
  )

rm(Africa)


# 5. Prepare iCOGS GWAS data
# ------------------------------------------------------------

iCOGS_miami <- iCOGS %>%
  transmute(
    Chr = as.numeric(chromosome),
    BP = as.numeric(base_pair_location),
    logp = -log10(p_value)
  ) %>%
  filter(
    !is.na(Chr),
    !is.na(BP),
    Chr >= 1,
    Chr <= 22,
    BP >= 0,
    is.finite(logp)
  )

rm(iCOGS)


# 6. Calculate shared chromosome offsets
# ------------------------------------------------------------
#
# Use the full GWAS datasets to calculate chromosome lengths
# so that both plots have consistent genomic coordinates.

chr_table <- bind_rows(
  Africa_miami %>% select(Chr, BP),
  iCOGS_miami %>% select(Chr, BP)
) %>%
  group_by(Chr) %>%
  summarise(
    chrLength = max(BP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Chr) %>%
  mutate(
    offset = c(0, head(cumsum(chrLength), -1)),
    chrEnd = offset + chrLength,
    chrMid = offset + chrLength / 2
  )


# 7. Define shared x-axis limits and chromosome labels
# ------------------------------------------------------------

total_length <- max(chr_table$chrEnd)

left_padding <- 0
right_padding <- total_length * 0.02

x_limits <- c(
  -left_padding,
  total_length + right_padding
)

x_breaks <- chr_table$chrMid

x_labels <- ifelse(
  chr_table$Chr %% 2 == 1 & chr_table$Chr != 19,
  as.character(chr_table$Chr),
  ""
)


# 8. Add cumulative genomic positions
# ------------------------------------------------------------

add_true_pos <- function(df, chr_table) {

  df <- df %>%
    mutate(
      Chr = as.numeric(as.character(Chr)),
      BP = as.numeric(as.character(BP))
    ) %>%
    select(-any_of(c("offset", "truePos")))

  chr_table_clean <- chr_table %>%
    transmute(
      Chr = as.numeric(as.character(Chr)),
      offset = as.numeric(as.character(offset))
    )

  df %>%
    left_join(chr_table_clean, by = "Chr") %>%
    mutate(truePos = BP + offset)
}


# Add genomic positions to Manhattan plot data
plotData <- add_true_pos(plotData, chr_table)

# Add genomic positions to Miami plot data
Africa_miami <- add_true_pos(Africa_miami, chr_table)
iCOGS_miami <- add_true_pos(iCOGS_miami, chr_table)


# 9. Generate Manhattan plot
# ------------------------------------------------------------

manPlot <- ggplot(
  plotData,
  aes(
    x = truePos,
    y = -log10(cum_avg_lfdr)
  )
) +

  geom_point(
    aes(color = "csmGmm replicated SNPs"),
    shape = 16
  ) +

  scale_color_manual(
    name = NULL,
    values = c(
      "csmGmm replicated SNPs" = "#009E73"
    ),
    breaks = "csmGmm replicated SNPs"
  ) +

  scale_x_continuous(
    name = "Chromosome",
    breaks = x_breaks,
    labels = x_labels,
    limits = x_limits,
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    name = "-log10(lfdr)"
  ) +

  theme_bw() +

  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )


# 10. Generate Miami plot
# ------------------------------------------------------------

miami_plot <- ggplot() +

  # Top: European ancestry (iCOGS)
  geom_point(
    data = iCOGS_miami,
    aes(
      x = truePos,
      y = logp,
      color = ifelse(
        Chr %% 2 == 0,
        "Even chromosome",
        "iCOGS"
      )
    ),
    alpha = 0.7,
    size = 0.8
  ) +

  # Bottom: African females
  geom_point(
    data = Africa_miami,
    aes(
      x = truePos,
      y = -logp,
      color = ifelse(
        Chr %% 2 == 0,
        "Even chromosome",
        "African"
      )
    ),
    alpha = 0.7,
    size = 0.8
  ) +

  scale_color_manual(
    name = NULL,
    values = c(
      "African" = "red",
      "iCOGS" = "blue",
      "Even chromosome" = "grey70",
      "csmGmm replicated SNPs" = "#009E73"
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    name = "Chromosome",
    breaks = x_breaks,
    labels = x_labels,
    limits = x_limits,
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    name = "-log10(p-value)",
    breaks = seq(-30, 80, 10)
  ) +

  coord_cartesian(
    ylim = c(-30, 80)
  ) +

  theme_cowplot() +

  theme(
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title = element_text(size = 24),
    plot.margin = margin(0, 5.5, 5.5, 5.5),
    legend.position = "none"
  )


# 11. Remove individual plot legends
# ------------------------------------------------------------

manPlot_noleg <- manPlot +
  theme(legend.position = "none")

miami_noleg <- miami_plot +
  theme(legend.position = "none")


# 12. Stack Manhattan and Miami plots
# ------------------------------------------------------------

main_plots <- plot_grid(
  manPlot_noleg,
  miami_noleg,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1.3)
)


# 13. Create shared legend
# ------------------------------------------------------------

legend_plot <- ggplot(
  data.frame(
    x = 1:3,
    y = 1:3,
    dataset = factor(
      c(
        "African",
        "iCOGS",
        "csmGmm replicated SNPs"
      ),
      levels = c(
        "African",
        "iCOGS",
        "csmGmm replicated SNPs"
      )
    )
  ),
  aes(x, y, color = dataset)
) +

  geom_point(size = 3) +

  scale_color_manual(
    name = NULL,
    values = c(
      "African" = "red",
      "iCOGS" = "blue",
      "csmGmm replicated SNPs" = "#009E73"
    ),
    breaks = c(
      "African",
      "iCOGS",
      "csmGmm replicated SNPs"
    )
  ) +

  theme_void() +

  theme(
    legend.position = "right",
    legend.text = element_text(size = 23)
  )


legend <- get_legend(legend_plot)


# 14. Combine plots and legend
# ------------------------------------------------------------

combined_plot <- plot_grid(
  main_plots,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.18)
)


# 15. Display combined figure
# ------------------------------------------------------------

print(combined_plot)


# 16. Save combined figure
# ------------------------------------------------------------

if (!dir.exists("results")) {
  dir.create("results")
}

ggsave(
  filename = "results/manhattan_miami_combined_new.png",
  plot = combined_plot,
  width = 27,
  height = 8,
  dpi = 300
)


# ============================================================
# End of example
# ============================================================
