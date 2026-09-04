## ---------------------------------------------------------------------------
## Figure 3 (SHAP interpretability), panels A-D, drawn in R.
##
##   A  mean absolute SHAP per feature, averaged over the 13 models
##   B  what each feature does to the prediction (shift plot, replaces the
##      beeswarm of the previous version)
##   C  contribution by feature family (C1) and by drug (C2)
##   D  co-variation among the leading features
##
## Heliyon formatting: Arial throughout, width kept under 20 cm, vector PDF plus
## 600 dpi uncompressed TIFF, no panel titles and no A/B/C letters.
##
## No panel titles and no A/B/C letters - those are added in BioRender.
## All text and numbers are bold.
##
## Panel B is no longer a beeswarm. With 20 test samples a beeswarm is mostly
## empty space and a heatmap of 20 columns is hard to read. It is now a shift
## plot: for each feature, the mean SHAP value when the feature is present (or
## high) against the mean when it is absent (or low), joined by an arrow. One
## glance gives the clinically meaningful quantity - how much carrying that
## feature moves the predicted risk, and in which direction.
##
## Input files, all in the working directory:
##   Figure_Data_All.xlsx        sheets Fig3A_Global_SHAP_Top15,
##                               Fig3C_Category_SHAP, Fig3C_Drug_SHAP,
##                               Fig3D_SHAP_Correlation
##   shap_gbm_test_full.csv      per-sample GBM SHAP values, all features
##   shap_gbm_test_features.csv  the value of each feature in the same 20
##                               samples, in the same order
##   shap_gbm_test_labels.csv    observed outcome and predicted probability
##                               for the same 20 samples
##
## Requires: ggplot2, dplyr, tidyr, readxl, ragg
## ---------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(ragg)

XLSX    <- "Figure_Data_All.xlsx"
B_SHAP  <- "shap_gbm_test_full.csv"
B_FEAT  <- "shap_gbm_test_features.csv"
B_LABEL <- "shap_gbm_test_labels.csv"
OUT     <- "."

POS <- "#C0392B"   # positive SHAP -> non-SVR
NEG <- "#2E8B7A"   # negative SHAP -> SVR
TXT <- "#1A1A1A"

FONT <- "Arial"    # Heliyon allows Arial, Helvetica, Courier, Times, Times New Roman, Symbol
BASE <- 11         # base font size; raise this if the panels will be shrunk
MAXW <- 7.8        # inches; 20 cm is the Heliyon maximum width

CAT_COL <- c("Genotype"            = "#2E7D8F",
             "Mutation count"      = "#A8C97F",
             "Drug RAS"            = "#E1873C",
             "Sequence similarity" = "#E6A0C4")

DIR_COL <- c("non-SVR" = POS, "SVR" = NEG, "Mixed" = "#A9A9A9")

CLASS_COL <- c("NS5A inhibitor"  = "#2E7D8F",
               "NS3/4A PI"       = "#E1873C",
               "NS5B polymerase" = "#A8C97F")

theme_panel <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      text               = element_text(family = FONT, face = "bold", colour = TXT),
      axis.text          = element_text(family = FONT, face = "bold", colour = TXT, size = base),
      axis.title         = element_text(family = FONT, face = "bold", colour = TXT, size = base + 1),
      plot.title         = element_blank(),
      axis.line          = element_line(colour = "#333333", linewidth = 0.6),
      axis.ticks         = element_line(colour = "#333333", linewidth = 0.6),
      panel.grid.major.x = element_line(colour = "#E0E0E0", linewidth = 0.35),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(family = FONT, face = "bold", size = base - 0.5),
      legend.key.size    = unit(0.42, "cm"),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.margin        = margin(10, 14, 8, 8)
    )
}

save_panel <- function(plot, name, w, h) {
  w <- min(w, MAXW)
  ## vector PDF for the journal
  ggsave(file.path(OUT, paste0(name, ".pdf")), plot,
         width = w, height = h, units = "in", device = cairo_pdf)
  ## 600 dpi PNG - BioRender accepts PNG but not TIFF, so use this for assembly
  agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
          units = "in", res = 600, background = "white")
  print(plot); dev.off()
  ## 600 dpi uncompressed TIFF for journal submission
  agg_tiff(file.path(OUT, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = 600, background = "white", compression = "none")
  print(plot); dev.off()
  cat("wrote", name, sprintf("(%.2f x %.2f in) pdf + png + tiff\n", w, h))
}


## Genotype_6r -> "Genotype 6r"; all_mutations_ns3 -> "All mutations NS3"
pretty_feature <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(ns3|ns5a|ns5b|NS3|NS5A|NS5B)\\b", "\\U\\1", x, perl = TRUE)
  x <- gsub("\\bmuts\\b", "RAS", x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

## ===========================================================================
## Panel A - global mean |SHAP| across the 13 models, with the number of
## models that rank each feature in their own top five
## ===========================================================================
A <- read_excel(XLSX, sheet = "Fig3A_Global_SHAP_Top15") %>%
  rename(feature = 1) %>%
  select(feature, mean_shap = Global_Mean_SHAP, top5 = Top5_Count, category = Category) %>%
  arrange(desc(mean_shap)) %>%
  mutate(lab = factor(pretty_feature(feature), levels = rev(pretty_feature(feature))))

xmaxA <- max(A$mean_shap) * 1.58

pA <- ggplot(A, aes(x = mean_shap, y = lab, fill = category)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.4f", mean_shap)),
            hjust = -0.15, size = 3.1, family = FONT, fontface = "bold", colour = "#333333") +
  geom_text(aes(x = xmaxA * 0.90, label = paste0(top5, "/13")),
            hjust = 0.5, size = 3.3, family = FONT, fontface = "bold", colour = TXT) +
  annotate("text", x = xmaxA * 0.88, y = nrow(A) + 1.15, size = 2.9, colour = "#555555",
           family = FONT, fontface = "bold", label = "models ranking\nthe feature\nin their top 5") +
  scale_x_continuous(limits = c(0, xmaxA), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = CAT_COL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_cartesian(clip = "off") +
  labs(x = "Global mean |SHAP value|", y = NULL) +
  theme_panel() +
  theme(plot.margin = margin(38, 20, 8, 8))

save_panel(pA, "Figure3A_importance", w = 7.0, h = 5.4)

## ===========================================================================
## Panel B - what each feature does to the prediction.
##
## For every feature the 20 test samples are split in two and the mean SHAP
## value of each group is plotted, joined by an arrow. The arrow therefore
## reads directly as "carrying this feature moves the predicted risk by this
## much, in this direction". Rows keep panel A's order.
##
## Split rule: value present versus absent when both groups have at least three
## samples (genotype dummies, mutation and RAS counts), otherwise a median split
## (the sequence-similarity features, which are continuous and never zero).
## ===========================================================================
Braw <- read.csv(B_SHAP, check.names = FALSE)
Vraw <- read.csv(B_FEAT, check.names = FALSE)
Lab  <- read.csv(B_LABEL, stringsAsFactors = FALSE)

id_col <- intersect(names(Braw), c("Sample", "sample_id", "sample"))[1]
want   <- tolower(A$feature)
gap    <- setdiff(want, tolower(setdiff(names(Braw), id_col)))
if (length(gap)) warning("panel B is missing SHAP for: ", paste(gap, collapse = ", "))

split_feature <- function(v) {
  if (sum(v > 0) >= 3 && sum(v == 0) >= 3) {
    list(g = v > 0, type = "present vs absent")
  } else {
    list(g = v > stats::median(v), type = "high vs low")
  }
}

rows <- lapply(seq_len(nrow(A)), function(i) {
  f  <- A$feature[i]
  sc <- names(Braw)[tolower(names(Braw)) == tolower(f)][1]
  vc <- names(Vraw)[tolower(names(Vraw)) == tolower(f)][1]
  if (is.na(sc) || is.na(vc)) return(NULL)
  sh <- Braw[[sc]]; v <- Vraw[[vc]]
  sp <- split_feature(v)
  data.frame(feature = f,
             lab     = as.character(A$lab[i]),
             type    = sp$type,
             hi_mean = mean(sh[sp$g]),  hi_n = sum(sp$g),
             lo_mean = mean(sh[!sp$g]), lo_n = sum(!sp$g))
})
B <- do.call(rbind, rows) %>%
  mutate(delta = hi_mean - lo_mean,
         dir   = ifelse(delta > 0, "Towards non-SVR", "Towards SVR"),
         lab   = factor(lab, levels = levels(A$lab)))

rng  <- range(c(B$hi_mean, B$lo_mean))
padL <- rng[1] - 0.10 * diff(rng)
padR <- rng[2] + 0.34 * diff(rng)

pB <- ggplot(B, aes(y = lab)) +
  geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.6) +
  geom_segment(aes(x = lo_mean, xend = hi_mean, yend = lab, colour = dir),
               linewidth = 1.5, lineend = "round",
               arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  geom_point(aes(x = lo_mean), size = 3.0, shape = 21, fill = "white",
             colour = "#555555", stroke = 1.1) +
  geom_point(aes(x = hi_mean, colour = dir), size = 3.4) +
  geom_text(aes(x = padR, label = sprintf("%+.3f", delta), colour = dir),
            hjust = 1, size = 3.1, family = FONT, fontface = "bold", show.legend = FALSE) +
  geom_text(aes(x = padR, y = nrow(B) + 1.0),
            label = "shift in mean\nSHAP value", hjust = 1, size = 2.9,
            family = FONT, fontface = "bold", colour = "#555555", check_overlap = TRUE) +
  scale_colour_manual(values = c("Towards non-SVR" = POS, "Towards SVR" = NEG),
                      name = NULL) +
  scale_x_continuous(limits = c(padL, padR), expand = expansion(mult = c(0.02, 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "Mean SHAP value (open circle: feature low or absent)", y = NULL) +
  theme_panel() +
  theme(plot.margin = margin(34, 18, 8, 8))

save_panel(pB, "Figure3B_shift", w = 7.0, h = 5.4)

## same plot without the feature names, for placing straight to the right of
## panel A so the two share one column of labels
pB_bare <- pB + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_panel(pB_bare, "Figure3B_shift_no_ylabels", w = 5.6, h = 5.4)

## ===========================================================================
## Panel C1 - contribution by feature family, coloured by net direction
## ===========================================================================
C1 <- read_excel(XLSX, sheet = "Fig3C_Category_SHAP") %>%
  rename(family = Category, total = Total_abs_SHAP, direction = Direction,
         nfeat = N_features) %>%
  arrange(total) %>%
  mutate(family = gsub("_", " ", family),
         family = factor(family, levels = family),
         direction = factor(direction, levels = c("non-SVR", "SVR", "Mixed")))

pC1 <- ggplot(C1, aes(x = total, y = family, fill = direction)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = paste0(nfeat, " features")),
            hjust = -0.12, size = 3.1, family = FONT, fontface = "bold", colour = "#444444") +
  scale_x_continuous(limits = c(0, max(C1$total) * 1.42),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = DIR_COL, name = NULL,
                    labels = c("Net effect towards non-SVR",
                               "Net effect towards SVR", "Mixed")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Total |SHAP| within family", y = NULL) +
  theme_panel() +
  theme(plot.margin = margin(10, 26, 8, 8))

save_panel(pC1, "Figure3C1_by_family", w = 6.8, h = 3.9)

## ===========================================================================
## Panel C2 - contribution by drug, coloured by drug class
## ===========================================================================
C2 <- read_excel(XLSX, sheet = "Fig3C_Drug_SHAP") %>%
  rename(drug = Drug, class = Drug_Class, shap = Mean_abs_SHAP) %>%
  arrange(desc(shap)) %>%
  mutate(drug = factor(drug, levels = drug))

pC2 <- ggplot(C2, aes(x = drug, y = shap, fill = class)) +
  geom_col(width = 0.72) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = CLASS_COL) +
  labs(x = NULL, y = "Mean |SHAP|") +
  theme_panel() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#E0E0E0", linewidth = 0.35),
        axis.text.x = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 1))

save_panel(pC2, "Figure3C2_by_drug", w = 6.4, h = 3.0)

## ===========================================================================
## Panel D - correlation between the SHAP vectors of the ten leading features
## ===========================================================================
Draw <- read_excel(XLSX, sheet = "Fig3D_SHAP_Correlation") %>% rename(feature_x = 1)
lev <- Draw$feature_x

D <- Draw %>%
  pivot_longer(-feature_x, names_to = "feature_y", values_to = "value") %>%
  mutate(value = ifelse(feature_x == feature_y, NA_real_, value),
         fx = factor(pretty_feature(feature_x), levels = pretty_feature(lev)),
         fy = factor(pretty_feature(feature_y), levels = rev(pretty_feature(lev))))

limD <- max(abs(D$value), na.rm = TRUE)

pD <- ggplot(D, aes(x = fx, y = fy, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value)),
                colour = abs(value) > 0.55 * limD),
            size = 3.0, family = FONT, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = TXT), na.value = TXT) +
  scale_fill_gradient2(low = NEG, mid = "#F7F7F7", high = "#B03A6E",
                       midpoint = 0, na.value = "#E8E8E8",
                       limits = c(-limD, limD), name = "Correlation of\nSHAP values") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_panel() +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = BASE - 1))

save_panel(pD, "Figure3D_covariation", w = 6.4, h = 5.4)

## ===========================================================================
## Verification
## ===========================================================================
cat("\n--- Panel A, top five ---\n")
print(as.data.frame(A[1:5, c("feature", "mean_shap", "top5", "category")]))
cat("\n--- Panel B ---\n")
print(as.data.frame(B[, c("feature", "type", "lo_n", "lo_mean", "hi_n", "hi_mean", "delta")]))
if (length(gap)) cat("MISSING:", paste(gap, collapse = ", "), "\n")
cat("row order identical in A and B:", identical(levels(A$lab), levels(B$lab)), "\n")
cat("observed outcomes in the 20 samples:", sum(Lab$True_Label == "Failure"),
    "non-SVR of", nrow(Lab), "\n")
cat("\n--- Panel C1 ---\n")
print(as.data.frame(C1[, c("family", "nfeat", "total", "direction")]))
cat("\n--- Panel D range ---\n")
print(range(D$value, na.rm = TRUE))
