# ============================================================================
# Figure 2 - STRATIFIED (OLS, with Dx): chord diagrams (AD | FTD | CN), EAA only
#   Input: aggregate results table (group-level statistics; NO participant data).
#   Reads sheets S_AD_complete / S_FTD_complete / S_CN_complete of EAA_primary_OLS_Dx.xlsx.
#   Outputs: BChord_AD_std.{png,svg}, BChord_FTD_std.{png,svg}, BChord_CN_std.{png,svg}.
# ============================================================================
if (!require("readxl"))   install.packages("readxl")
if (!require("circlize")) install.packages("circlize")
if (!require("dplyr"))    install.packages("dplyr")
if (!require("stringr"))  install.packages("stringr")
if (!require("svglite"))  install.packages("svglite")
library(readxl); library(circlize); library(dplyr); library(stringr); library(svglite); library(grid)

# ---- Paths (portable; no personal paths) ----
# Run this script from the folder that holds EAA_primary_OLS_Dx.xlsx (in RStudio:
# Session > Set Working Directory > To Source File Location). Figures are written here.
WORKDIR   <- getwd()
if (!dir.exists(WORKDIR)) dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
EAA_FILE  <- file.path(WORKDIR, "EAA_primary_OLS_Dx.xlsx")   # S_AD_complete / S_FTD_complete / S_CN_complete
FONT <- "Helvetica"
nonsig_color   <- "#C2C2C2"
C_FDR <- "#D62728"; C_NOM <- "#E8887C"; C_NS <- "#C2C2C2"
nonsig_ribbon  <- "#F2F2F2"
family_headers <- c("First/Second Generation Clocks", "Next Generation Clocks", "Organs/Systems Clocks")
fam_levels     <- c("First/Second Generation", "Next Generation", "Organs/Systems")
family_base    <- c("First/Second Generation"="#9C27B0","Next Generation"="#F57B00","Organs/Systems"="#009E73")
fam_label      <- c("First/Second Generation"="First/Second","Next Generation"="Next","Organs/Systems"="Organs/Systems")

shades <- function(base, n, lite=0.55, dark=0.40) {
  if (n == 1) return(base)
  brgb <- col2rgb(base)[,1] / 255
  mixw <- function(a) grDevices::rgb(brgb[1]+(1-brgb[1])*a, brgb[2]+(1-brgb[2])*a, brgb[3]+(1-brgb[3])*a)
  mixb <- function(a) grDevices::rgb(brgb[1]*(1-a), brgb[2]*(1-a), brgb[3]*(1-a))
  ts <- seq(lite, -dark, length.out = n)
  sapply(ts, function(t) if (t >= 0) mixw(t) else mixb(-t))
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
           sig_p     = P_FDR_num,
           Label     = str_replace(Measure, "^EAA ", ""))
  d$Family <- factor(d$Family, levels = fam_levels)
  d <- d %>% filter(!is.na(beta_std), !is.na(Label), Label != "") %>% arrange(Family)
  d$ClockColor <- NA_character_
  for (fam in levels(d$Family)) {
    idx <- which(d$Family == fam); if (!length(idx)) next
    d$ClockColor[idx] <- shades(family_base[[fam]], length(idx))
  }
  d$is_sig    <- (!is.na(d$sig_p) & d$sig_p < 0.05) | (!is.na(d$P_num) & d$P_num < 0.05)
  d$LinkColor <- ifelse(d$is_sig, d$ClockColor, nonsig_ribbon)
  d
}

dAD  <- prep_sheet_std(EAA_FILE, "S_AD_complete")
dFTD <- prep_sheet_std(EAA_FILE, "S_FTD_complete")
dCN  <- prep_sheet_std(EAA_FILE, "S_CN_complete")

# ---- Chord builder ----
make_chord <- function(d, fam_text = TRUE) {
  circos.clear()
  d <- d[is.finite(d$beta_std) & !is.na(d$Label) & d$Label != "" & !is.na(d$ClockColor), ]
  d <- d[!duplicated(d$Label), ]
  vals <- pmax(abs(d$beta_std), 0.01)
  mat <- matrix(vals, nrow = 1, dimnames = list("NPS", d$Label))
  grid.col <- c("NPS" = "#F0E442", setNames(d$ClockColor, d$Label))
  link.col <- d$LinkColor
  gaps <- setNames(rep(1.2, ncol(mat) + 1), c(colnames(mat), "NPS")); gaps["NPS"] <- 8.5
  circos.par(gap.after = gaps, start.degree = 17.5, points.overflow.warning = FALSE,
             canvas.xlim = c(-1.05, 1.05), canvas.ylim = c(-1.05, 1.05),
             circle.margin = c(0.05, 0.05, 0.05, 0.05))
  chordDiagram(mat, grid.col = grid.col, col = link.col, directional = 0,
               transparency = 0, annotationTrack = "grid", link.border = NA,
               reduce = -1, preAllocateTracks = list(track.height = 0.10))
  circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    sn <- get.cell.meta.data("sector.index")
    xl <- get.cell.meta.data("xlim"); yl <- get.cell.meta.data("ylim")
    if (sn == "NPS") circos.text(mean(xl), yl[2] + 0.05, "NPS",
                                 facing = "inside", niceFacing = TRUE, cex = 0.58, font = 1)
  })
  
  if (fam_text) {
    for (fam in levels(d$Family)) {
      secs <- as.character(d$Label[d$Family == fam]); if (!length(secs)) next
      highlight.sector(secs, track.index = 1, col = NA, text = fam_label[[fam]],
                       cex = 0.50, text.col = "black", facing = "bending.outside", niceFacing = TRUE)
    }
  }
  circos.clear()
}

# ---- Export ----
H_CHORD <- 90
close_all_dev <- function() { while (!is.null(dev.list())) try(dev.off(), silent = TRUE) }
save_dev <- function(fname, w_mm, h_mm, draw_fun, svg = FALSE) {
  close_all_dev()
  path <- file.path(getwd(), fname)
  if (file.exists(path)) try(file.remove(path), silent = TRUE)
  if (svg) svglite::svglite(filename = path, width = w_mm/25.4, height = h_mm/25.4, pointsize = 12)
  else     grDevices::png(filename = path, width = w_mm/25.4, height = h_mm/25.4, units = "in", res = 450, pointsize = 12)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  circos.clear(); draw_fun(); circos.clear()
}
for (s in list(list("AD", dAD), list("FTD", dFTD), list("CN", dCN))) {
  tag <- s[[1]]; d <- s[[2]]
  save_dev(sprintf("BChord_%s_std.png", tag), 70, H_CHORD, function() make_chord(d, TRUE), svg = FALSE)
  save_dev(sprintf("BChord_%s_std.svg", tag), 70, H_CHORD, function() make_chord(d, TRUE), svg = TRUE)
}
close_all_dev()
cat("\nStratified chords (AD, FTD, CN), with Dx, EAA only, saved in:", getwd(), "\n")
