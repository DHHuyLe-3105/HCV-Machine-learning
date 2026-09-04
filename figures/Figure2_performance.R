## ---------------------------------------------------------------------------
## Figure 2 - model performance across the five feature sets.
##   A  AUC heatmap, 13 classifiers x 5 feature sets
##   B  whole-genome versus three-gene predicted probabilities, per model
##   C  best AUC achieved by each feature set
##   D  ROC curves of the three best models on the RAS feature set
##
## Heliyon formatting: Arial throughout, single column, width kept under 20 cm,
## vector PDF plus 600 dpi uncompressed TIFF. No panel titles and no A/B/C
## letters - add those in BioRender.
##
## Input: Figure_Data_All.xlsx
## Requires: ggplot2, dplyr, tidyr, readxl, ragg
## ---------------------------------------------------------------------------

library(ggplot2); library(dplyr); library(tidyr); library(readxl); library(ragg)
library(scales)

XLSX <- "Figure_Data_All.xlsx"; OUT <- "."
FONT <- "Arial"      # Heliyon allows Arial, Helvetica, Courier, Times, Times New Roman, Symbol
BASE <- 11
TXT  <- "#1A1A1A"
MAXW <- 7.8          # inches; 20 cm is the Heliyon maximum width

theme_panel <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      text             = element_text(family = FONT, face = "bold", colour = TXT),
      axis.text        = element_text(family = FONT, face = "bold", colour = TXT, size = base),
      axis.title       = element_text(family = FONT, face = "bold", colour = TXT, size = base + 1),
      plot.title       = element_blank(),
      axis.line        = element_line(colour = "#333333", linewidth = 0.6),
      axis.ticks       = element_line(colour = "#333333", linewidth = 0.6),
      legend.title     = element_text(family = FONT, face = "bold", size = base - 1),
      legend.text      = element_text(family = FONT, face = "bold", size = base - 1),
      strip.text       = element_text(family = FONT, face = "bold", size = base - 1),
      strip.background = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(10, 14, 8, 8)
    )
}

save_panel <- function(p, name, w, h) {
  w <- min(w, MAXW)
  ## vector PDF for the journal
  ggsave(file.path(OUT, paste0(name, ".pdf")), p,
         width = w, height = h, units = "in", device = cairo_pdf)
  ## 600 dpi PNG - BioRender accepts PNG but not TIFF, so use this for assembly
  agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
          units = "in", res = 600, background = "white")
  print(p); dev.off()
  ## 600 dpi uncompressed TIFF for journal submission
  agg_tiff(file.path(OUT, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = 600, background = "white", compression = "none")
  print(p); dev.off()
  cat("wrote", name, sprintf("(%.2f x %.2f in) pdf + png + tiff\n", w, h))
}


SETS <- c("Nu_WG", "Nu_3TG", "Aa_WG", "Aa_3TG", "RAS")

## Short codes on the axes so the labels cannot collide. They must be expanded
## in the figure caption, exactly as printed by the check at the end of this
## script: Nu, nucleotide; Aa, amino acid; WG, whole genome; 3TG, the three
## target genes NS3, NS5A and NS5B; RAS, the 95 resistance and subtype features.
SETLAB <- c(Nu_WG = "Nu_WG", Nu_3TG = "Nu_3TG", Aa_WG = "Aa_WG",
            Aa_3TG = "Aa_3TG", RAS = "RAS")
SETFULL <- c(Nu_WG = "nucleotide, whole genome", Nu_3TG = "nucleotide, three genes",
             Aa_WG = "amino acid, whole genome", Aa_3TG = "amino acid, three genes",
             RAS = "resistance and subtype features")

## ------------------------------------------------------------ A: heatmap
H <- read_excel(XLSX, sheet = "Fig2A_Heatmap_AUC") %>%
  pivot_longer(-Model, names_to = "set", values_to = "auc") %>%
  mutate(set = sub("_AUC$", "", set),
         set = factor(SETLAB[set], levels = SETLAB[SETS]))
ordM <- H %>% group_by(Model) %>% summarise(m = mean(auc)) %>% arrange(m)
H$Model <- factor(H$Model, levels = ordM$Model)

pA <- ggplot(H, aes(x = set, y = Model, fill = auc)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", auc), colour = auc > 0.90),
            family = FONT, fontface = "bold", size = 3.0, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = TXT)) +
  scale_fill_gradientn(
    colours = c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272",
                "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    values  = scales::rescale(c(0, 0.30, 0.55, 0.70, 0.80, 0.86, 0.91, 0.94, 1.00)),
    limits = c(0, 1), name = "AUC") +
  guides(fill = guide_colourbar(barwidth = unit(0.4, "cm"), barheight = unit(3.4, "cm"))) +
  labs(x = NULL, y = NULL) +
  theme_panel() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(family = FONT, face = "bold", size = BASE + 0.5))
save_panel(pA, "Figure2A_heatmap", 6.8, 5.0)

## ------------------------------- B: whole genome versus three genes
## keys must be the model names EXACTLY as they appear in the t-test sheet,
## otherwise the P values do not join; values are the sheet-name suffixes
sheets <- c("GBM"                 = "GBM",
            "Stacking Classifier" = "Stacking",
            "Voting Classifier"   = "Voting_C",
            "Random Forest"       = "Random_F",
            "Logistic Regression" = "Logistic")

## shorter text for the facet strips only
STRIP <- c("GBM" = "GBM", "Stacking Classifier" = "Stacking",
           "Voting Classifier" = "Voting", "Random Forest" = "Random forest",
           "Logistic Regression" = "Logistic regression")
grab <- function(prefix, enc) {
  bind_rows(lapply(names(sheets), function(lab) {
    s <- read_excel(XLSX, sheet = paste0("Fig2B_", prefix, "_", sheets[[lab]]))
    names(s)[3:4] <- c("WG", "TG")
    s %>% select(Sample, True_Label, WG, TG) %>%
      mutate(model = lab, encoding = enc)
  }))
}
## the t-test sheet writes "Amino Acid" with a capital A, so both tables are
## normalised before faceting - otherwise the facets do not match and only one
## P value is drawn
norm_enc <- function(x) ifelse(grepl("amino", tolower(x)), "Amino acid", "Nucleotide")

B <- bind_rows(grab("Am", "Amino acid"), grab("Nu", "Nucleotide")) %>%
  pivot_longer(c(WG, TG), names_to = "scope", values_to = "prob") %>%
  mutate(scope = factor(ifelse(scope == "WG", "Whole genome", "Three genes"),
                        levels = c("Whole genome", "Three genes")),
         encoding = factor(norm_enc(encoding), levels = c("Amino acid", "Nucleotide")),
         model = factor(model, levels = names(sheets)))

TT <- read_excel(XLSX, sheet = "Fig2B_WGvs3TG_ttest") %>%
  filter(Model %in% names(sheets)) %>%
  transmute(model = factor(Model, levels = names(sheets)),
            encoding = factor(norm_enc(Encoding), levels = c("Amino acid", "Nucleotide")),
            lab = ifelse(p_value < 0.001, "P < 0.001", sprintf("P = %.3f", p_value)))
if (nrow(TT) != 10) {
  stop("expected 10 P values, got ", nrow(TT),
       ". Model names in Fig2B_WGvs3TG_ttest are: ",
       paste(sort(unique(read_excel(XLSX, sheet = "Fig2B_WGvs3TG_ttest")$Model)),
             collapse = ", "),
       " - the names(sheets) keys must match these exactly.")
}

pB <- ggplot(B, aes(x = scope, y = prob, fill = scope)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.30, colour = TXT,
               linewidth = 0.5) +
  geom_jitter(aes(colour = scope), width = 0.15, size = 1.0, alpha = 0.75,
              show.legend = FALSE) +
  geom_text(data = TT, aes(x = 1.5, y = 108, label = lab), inherit.aes = FALSE,
            family = FONT, fontface = "bold", size = 2.8, colour = TXT) +
  facet_grid(encoding ~ model, labeller = labeller(model = STRIP)) +
  scale_fill_manual(values = c("Whole genome" = "#2E8B7A", "Three genes" = "#E07B39"),
                    name = NULL) +
  scale_colour_manual(values = c("Whole genome" = "#1F6F63", "Three genes" = "#B35A00")) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 115)) +
  labs(x = NULL, y = "Predicted probability of non-SVR (%)") +
  theme_panel() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = BASE - 1.5),
        legend.position = "none",
        panel.grid.major.y = element_line(colour = "#EDEDED", linewidth = 0.3))
save_panel(pB, "Figure2B_WG_vs_3TG", 7.8, 4.6)

## ------------------------------------------------- C: best AUC per feature set
C <- read_excel(XLSX, sheet = "Fig2C_BestAUC_BarChart") %>%
  mutate(lab = factor(SETLAB[Feature_Set], levels = SETLAB[SETS]),
         hi = Feature_Set == "RAS")

pC <- ggplot(C, aes(x = lab, y = Best_AUC, fill = hi)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f\n%s", Best_AUC, Best_Model)),
            vjust = -0.25, family = FONT, fontface = "bold", size = 3.0, colour = TXT) +
  geom_text(aes(y = 0.055, label = paste0(format(N_Features, big.mark = ","), "\nfeatures"),
                colour = hi), family = FONT, fontface = "bold", size = 2.8,
            show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#333333")) +
  scale_fill_manual(values = c(`TRUE` = "#A50F15", `FALSE` = "#FC9272")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Best AUC on the held-out split") +
  theme_panel() +
  theme(panel.grid.major.y = element_line(colour = "#EDEDED", linewidth = 0.3))
save_panel(pC, "Figure2C_best_AUC", 6.4, 4.2)

## ---------------------------------------------------------- D: ROC curves
roc <- function(sheet, label) {
  read_excel(XLSX, sheet = sheet) %>% select(FPR, TPR) %>%
    filter(!is.na(FPR), !is.na(TPR)) %>% mutate(model = label)
}
## Youden J operating point: the threshold maximising sensitivity + specificity - 1
youden <- function(df) df %>% group_by(model) %>%
  slice_max(TPR - FPR, n = 1, with_ties = FALSE) %>% ungroup()
R <- bind_rows(roc("Fig2D_ROC_Stacking_Cla", "Stacking classifier (0.945)"),
               roc("Fig2D_ROC_Voting_Class", "Voting classifier (0.945)"),
               roc("Fig2D_ROC_GBM",          "Gradient boosting (0.944)"))

pD <- ggplot(R, aes(x = FPR, y = TPR, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "#BBBBBB",
              linetype = "22", linewidth = 0.6) +
  geom_step(linewidth = 1.1, direction = "hv") +
  geom_point(data = youden(R), aes(x = FPR, y = TPR, colour = model),
             shape = 8, size = 3.2, stroke = 1.1, show.legend = FALSE) +
  scale_colour_manual(values = c("Stacking classifier (0.945)" = "#5B4B9E",
                                 "Voting classifier (0.945)"   = "#1F6F63",
                                 "Gradient boosting (0.944)"   = "#C0392B"), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(0.01)) +
  coord_fixed() +
  labs(x = "1 - specificity", y = "Sensitivity") +
  theme_panel() +
  theme(legend.position = c(0.63, 0.20),
        legend.background = element_rect(fill = "white", colour = "#DDDDDD"))
save_panel(pD, "Figure2D_ROC", 5.2, 5.2)

## ------------------------------------------------------------------ checks
cat("\nAxis codes to define in the Figure 2 caption:\n")
for (k in SETS) cat("  ", k, "=", SETFULL[[k]], "\n")
cat("\nBest model per feature set:\n")
print(as.data.frame(C[, c("Feature_Set", "N_Features", "Best_Model", "Best_AUC")]))
cat("\nModels with AUC below 0.5 (inverted or failed fits):\n")
print(as.data.frame(H[H$auc < 0.5, c("Model", "set", "auc")]))
cat("\nSignificant WG vs 3TG differences:\n")
print(as.data.frame(read_excel(XLSX, sheet = "Fig2B_WGvs3TG_ttest") %>%
                      filter(Significant_p005 == "Yes") %>%
                      select(Encoding, Model, p_value)))
