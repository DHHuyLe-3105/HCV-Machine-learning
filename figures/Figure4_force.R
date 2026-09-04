## ---------------------------------------------------------------------------
## Figure 4 - individual prediction decomposition (force plots) for four
## representative test cases, redrawn in R from the deposited SHAP values.
##
## Each panel is one case: features ordered by absolute contribution, bars
## coloured by direction, the summed SHAP value stated in the corner. The four
## panels are written separately so they can be arranged in BioRender.
##
## Heliyon formatting: Arial, width under 20 cm, vector PDF plus 600 dpi
## uncompressed TIFF, no panel titles and no A/B/C letters.
##
## Input: Figure_Data_All.xlsx, sheet Fig4_Force_Plots
## Requires: ggplot2, dplyr, readxl, ragg
## ---------------------------------------------------------------------------

library(ggplot2); library(dplyr); library(readxl); library(ragg)

XLSX <- "Figure_Data_All.xlsx"; OUT <- "."
FONT <- "Arial"
BASE <- 11
MAXW <- 7.8
TXT  <- "#1A1A1A"
POS  <- "#C0392B"   # pushes towards non-SVR
NEG  <- "#2E8B7A"   # pushes towards SVR

theme_panel <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      text               = element_text(family = FONT, face = "bold", colour = TXT),
      axis.text          = element_text(family = FONT, face = "bold", colour = TXT, size = base),
      axis.title         = element_text(family = FONT, face = "bold", colour = TXT, size = base + 1),
      plot.title         = element_blank(),
      axis.line          = element_line(colour = "#333333", linewidth = 0.6),
      axis.ticks         = element_line(colour = "#333333", linewidth = 0.6),
      panel.grid.major.x = element_line(colour = "#EDEDED", linewidth = 0.35),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(family = FONT, face = "bold", size = base - 1),
      legend.key.size    = unit(0.4, "cm"),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.margin        = margin(10, 14, 8, 8)
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


pretty_feature <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(ns3|ns5a|ns5b|NS3|NS5A|NS5B)\\b", "\\U\\1", x, perl = TRUE)
  x <- gsub("\\bmuts\\b", "RAS", x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

F <- read_excel(XLSX, sheet = "Fig4_Force_Plots") %>%
  mutate(lab = pretty_feature(Feature),
         dir = ifelse(SHAP_Value > 0, "Towards non-SVR", "Towards SVR"))

cases <- F %>% distinct(Sample_ID, Sample_Label, Total_SHAP) %>% arrange(desc(Total_SHAP))
lim <- max(abs(F$SHAP_Value)) * 1.45

for (k in seq_len(nrow(cases))) {
  cs <- cases[k, ]
  dat <- F %>% filter(Sample_ID == cs$Sample_ID) %>%
    arrange(Abs_SHAP) %>%
    mutate(lab = factor(lab, levels = lab))
  arrow_txt <- sprintf("Sum of SHAP = %+.3f  \u2192  %s",
                       cs$Total_SHAP, ifelse(cs$Total_SHAP > 0, "non-SVR", "SVR"))
  p <- ggplot(dat, aes(x = SHAP_Value, y = lab, fill = dir)) +
    geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.6) +
    geom_col(width = 0.66) +
    geom_text(aes(label = sprintf("%+.3f", SHAP_Value),
                  hjust = ifelse(SHAP_Value > 0, -0.18, 1.18)),
              family = FONT, fontface = "bold", size = 3.0, colour = "#333333") +
    annotate("label", x = -lim, y = nrow(dat) + 0.85, hjust = 0, label = arrow_txt,
             family = FONT, fontface = "bold", size = 3.0, colour = TXT,
             fill = ifelse(cs$Total_SHAP > 0, "#FDF2F0", "#EDF6F3"),
             label.size = 0.5) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.5))) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("Towards non-SVR" = POS, "Towards SVR" = NEG),
                      name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(-lim, lim)) +
    labs(x = "SHAP value", y = NULL) +
    theme_panel() +
    theme(plot.margin = margin(16, 14, 8, 8))
  save_panel(p, sprintf("Figure4%s_%s", LETTERS[k],
                        gsub("[^A-Za-z0-9]+", "_", cs$Sample_Label)), 6.4, 3.9)
}

cat("\ncases drawn:\n"); print(as.data.frame(cases))
cat("\nlargest single contribution per case:\n")
print(as.data.frame(F %>% group_by(Sample_Label) %>%
                      slice_max(Abs_SHAP, n = 1) %>%
                      select(Sample_Label, Feature, SHAP_Value)))
