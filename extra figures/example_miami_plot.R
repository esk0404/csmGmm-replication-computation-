# ============================================================
# Example: Miami Plot for Breast Cancer GWAS Replication
# ============================================================
#
# This example generates a Miami plot comparing GWAS
# association statistics from:
#
#   1. African females (GCST90551892)
#   2. European females (iCOGS, GCST010098)
#
# Both datasets use genome build 37.
#
# Data download instructions are provided in the data/ folder.
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
#
# Download the iCOGS dataset (GCST010098) following the
# instructions in the data/ folder.
#
# Place the downloaded file in the data/ directory.
#
# Expected files:
#   data/GCST010098.tsv
#   data/GCST90551892.tsv.gz


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


# 3. Prepare African GWAS data
# ------------------------------------------------------------

Africa_miami <- Africa %>%
  transmute(
    Chr = as.integer(chromosome),
    BP = base_pair_location,
    logp = -log10(p_value)
  ) %>%
  filter(
    !is.na(Chr),
    Chr >= 1,
    Chr <= 22,
    !is.na(BP),
    is.finite(logp)
  )

rm(Africa)


# 4. Prepare iCOGS GWAS data
# ------------------------------------------------------------

iCOGS_miami <- iCOGS %>%
  transmute(
    Chr = as.integer(chromosome),
    BP = base_pair_location,
    logp = -log10(p_value)
  ) %>%
  filter(
    !is.na(Chr),
    Chr >= 1,
    Chr <= 22,
    !is.na(BP),
    is.finite(logp)
  )

rm(iCOGS)


# 5. Calculate chromosome offsets
# ------------------------------------------------------------

chr_table <- bind_rows(
  Africa_miami %>% select(Chr, BP),
  iCOGS_miami %>% select(Chr, BP)
) %>%
  group_by(Chr) %>%
  summarise(
    chrLength = max(BP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Chr)

chrOffsets <- cumsum(chr_table$chrLength)
names(chrOffsets) <- chr_table$Chr


# 6. Calculate cumulative genomic positions
# ------------------------------------------------------------

add_true_pos <- function(df, chrOffsets) {
  df %>%
    arrange(Chr, BP) %>%
    mutate(
      truePos = BP + ifelse(
        Chr == 1,
        0,
        chrOffsets[as.character(Chr - 1)]
      )
    )
}

Africa_miami <- add_true_pos(Africa_miami, chrOffsets)
iCOGS_miami <- add_true_pos(iCOGS_miami, chrOffsets)


# 7. Define chromosome axis labels
# ------------------------------------------------------------

xBreaks <- chrOffsets[as.character(1:22)]

xBreaksLabs <- ifelse(
  (1:22) %% 2 == 0,
  "",
  as.character(1:22)
)


# 8. Generate Miami plot
# ------------------------------------------------------------

miami_plot <- ggplot() +

  # Top: European females (iCOGS)
  geom_point(
    data = iCOGS_miami,
    aes(x = truePos, y = logp),
    color = ifelse(
      iCOGS_miami$Chr %% 2 == 0,
      "grey70",
      "blue"
    ),
    alpha = 0.7,
    size = 0.8
  ) +

  # Bottom: African females
  geom_point(
    data = Africa_miami,
    aes(x = truePos, y = -logp),
    color = ifelse(
      Africa_miami$Chr %% 2 == 0,
      "grey70",
      "red"
    ),
    alpha = 0.7,
    size = 0.8
  ) +

  scale_x_continuous(
    name = "Chromosome",
    breaks = xBreaks,
    labels = xBreaksLabs
  ) +

  scale_y_continuous(
    name = "-log10(p-value)",
    limits = c(-30, 80),
    breaks = seq(-30, 80, 10)
  ) +

  theme_cowplot() +

  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )


# Display plot
print(miami_plot)


# 9. Save Miami plot
# ------------------------------------------------------------

if (!dir.exists("results")) {
  dir.create("results")
}

ggsave(
  filename = "results/miami_plot_rep_new.png",
  plot = miami_plot,
  width = 14,
  height = 6,
  dpi = 300
)
