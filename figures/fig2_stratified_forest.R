# ============================================================================
# Figure 2 - STRATIFIED (OLS, with Dx): forest plot per panel (AD | FTD | CN), EAA only
#   Input: aggregate results table (group-level statistics; NO participant data).
#   Reads sheets S_AD_complete / S_FTD_complete / S_CN_complete of EAA_primary_OLS_Dx.xlsx.
#   Outputs: BForest_AD_std / BForest_FTD_std / BForest_CN_std ({png,svg}).
# ============================================================================
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("stringr")) install.packages("stringr")
if (!require("svglite")) install.packages("svglite")
library(readxl); library(dplyr); library(ggplot2); library(stringr); library(svglite); library(grid)

# ---- Paths (portable; no personal paths) ----
# Run this script from the folder that holds EAA_primary_OLS_Dx.xlsx (in RStudio:
# Session > Set Working Directory > To Source File Location). Figures are written here.
WORKDIR <- getwd()
while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
if (!dir.exists(WORKDIR)) dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

EAA_FILE <- file.path(WORKDIR, "EAA_primary_OLS_Dx.xlsx")   # S_AD_complete / S_FTD_complete / S_CN_complete
FONT <- "Helvetica"
C_FDR <- "#D62728"; C_NOM <- "#E8887C"; nonsig_color <- "#C2C2C2"
family_headers <- c("First/Second Generation Clocks","Next Generation Clocks","Organs/Systems Clocks")
fam_levels     <- c("First/Second Generation","Next Generation","Organs/Systems")
family_base    <- c("First/Second Generation"="#9C27B0","Next Generation"="#F57B00","Organs/Systems"="#009E73")

shades <- function(base, n, lite=0.55, dark=0.40) {
  if (n == 1) return(base)
  brgb <- col2rgb(base)[,1] / 255
  mixw <- function(a) grDevices::rgb(brgb[1]+(1-brgb[1])*a, brgb[2]+(1-brgb[2])*a, brgb[3]+(1-brgb[3])*a)
  mixb <- function(a) grDevices::rgb(brgb[1]*(1-a), brgb[2]*(1-a), brgb[3]*(1-a))
  ts <- seq(lite, -dark, length.out = n)
  sapply(ts, function(t) if (t >= 0) mixw(t) else mixb(-t))
}
sig_color <- function(pfdr, pnom) {
  if (!is.na(pfdr) && pfdr < 0.05) return(C_FDR)
  if (!is.na(pnom) && pnom < 0.05) return(C_NOM)
  return(nonsig_color)
}

# ---- EAA stratified: Measure | beta | 95% CI (std) | P | P-FDR | Cohen's d | Partial R2 ----
prep_sheet_std <- function(file_path, sheet_name) {
  raw <- read_excel(file_path, sheet = sheet_name)
  raw <- raw[, 1:7]
  names(raw) <- c("Measure","beta_std","CI","P","P_FDR","Cohen_d","R2")
  cf <- NA; fv <- character(nrow(raw))
  for (i in seq_len(nrow(raw))) {
    val <- raw$Measure[i]
    if (!is.na(val) && val %in% family_headers) {
      cf <- dplyr::case_when(
        val == "First/Second Generation Clocks" ~ "First/Second Generation",
        val == "Next Generation Clocks"         ~ "Next Generation",
        val == "Organs/Systems Clocks"          ~ "Organs/Systems")
    }
    fv[i] <- cf
  }
  raw$Family <- fv
  d <- raw %>%
    filter(!is.na(Measure), !(Measure %in% family_headers)) %>%
    mutate(beta_std = as.numeric(beta_std),
           P_num     = as.numeric(str_extract(as.character(P), "^[0-9.]+")),
           P_FDR_num = suppressWarnings(as.numeric(str_extract(as.character(P_FDR),"^[0-9.]+"))),
           CI_L_std  = as.numeric(str_extract(CI, "(?<=\\[)[-0-9.]+")),
           CI_U_std  = as.numeric(str_extract(CI, "(?<=,\\s)[-0-9.]+")),
           sig_p     = P_FDR_num,
           Label     = str_replace(Measure, "^EAA ", ""),
           Label     = ifelse(Label %in% c("IntrinsicCapacity","DunedinPACE"), paste0("* ", Label), Label)) %>%
    mutate(.l = pmin(CI_L_std, CI_U_std), .u = pmax(CI_L_std, CI_U_std), CI_L_std = .l, CI_U_std = .u) %>% select(-.l,-.u)
  d$Family <- factor(d$Family, levels = fam_levels)
  d <- d %>% arrange(Family)
  d$ClockColor <- NA_character_
  for (fam in levels(d$Family)) {
    idx <- which(d$Family == fam); if (!length(idx)) next
    d$ClockColor[idx] <- shades(family_base[[fam]], length(idx))
  }
  d$FillColA  <- mapply(sig_color, d$sig_p, d$P_num)
  d$FillSolid <- d$FillColA
  d
}

dAD  <- prep_sheet_std(EAA_FILE, "S_AD_complete")
dFTD <- prep_sheet_std(EAA_FILE, "S_FTD_complete")
dCN  <- prep_sheet_std(EAA_FILE, "S_CN_complete")

# Forest builder (per-panel scales/breaks)
make_forest_custom <- function(d, x_lo, x_hi, brks, labels = TRUE) {
  df <- d %>% mutate(Label = factor(Label, levels = rev(Label)))
  span <- x_hi - x_lo; sq_w <- span * 0.048; sq_x <- x_lo - (span * 0.040)
  p <- ggplot(df, aes(y = Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey50") +
    geom_errorbar(aes(xmin = CI_L_std, xmax = CI_U_std, colour = FillColA), orientation = "y", width = 0.28, linewidth = 0.7) +
    geom_point(aes(x = beta_std, colour = FillColA), size = 2.0) +
    scale_colour_identity()
  if (labels) p <- p + geom_tile(aes(x = sq_x, fill = ClockColor), width = sq_w, height = 0.68, colour = "black", linewidth = 0.1)
  p + scale_fill_identity() +
    scale_x_continuous(expression(paste("Standardized Effect (", beta[std], ")")), breaks = brks) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") + labs(y = NULL) +
    theme_classic(base_size = 7) +
    theme(text = element_text(family = FONT),
          axis.text.y  = if (labels) element_text(size = 5.5, colour = "black", margin = margin(r = 4)) else element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x  = element_text(size = 4.5, colour = "black"), axis.title.x = element_text(size = 5),
          axis.line.x  = element_line(linewidth = 0.15), axis.ticks.x = element_line(linewidth = 0.15),
          plot.margin  = margin(2, 2, 2, if (labels) 20 else 2))
}

fAD  <- make_forest_custom(dAD,  -0.25, 0.40, seq(-0.2, 0.4, 0.2),  labels = TRUE)
fCN  <- make_forest_custom(dCN,  -0.30, 0.40, seq(-0.2, 0.4, 0.2),  labels = FALSE)
fFTD <- make_forest_custom(dFTD, -0.60, 0.60, seq(-0.5, 0.5, 0.25), labels = FALSE)

H_FOREST <- 62
save_gg <- function(fname, plot_obj, w_mm, h_mm, svg=FALSE) {
  while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
  p <- file.path(getwd(), fname); if (file.exists(p)) try(file.remove(p), silent=TRUE)
  if (svg) ggsave(p, plot_obj, width=w_mm/25.4, height=h_mm/25.4, device=svglite::svglite, bg="white")
  else     ggsave(p, plot_obj, width=w_mm/25.4, height=h_mm/25.4, dpi=450, bg="white")
  while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
}
save_gg("BForest_AD_std.png",  fAD,  55, H_FOREST)
save_gg("BForest_CN_std.png",  fCN,  35, H_FOREST)
save_gg("BForest_FTD_std.png", fFTD, 35, H_FOREST)
save_gg("BForest_AD_std.svg",  fAD,  55, H_FOREST, svg=TRUE)
save_gg("BForest_CN_std.svg",  fCN,  35, H_FOREST, svg=TRUE)
save_gg("BForest_FTD_std.svg", fFTD, 35, H_FOREST, svg=TRUE)
cat("Saved stratified forests (with Dx, EAA only) in:", getwd(), "\n")
