# ============================================================================
# Figure 2 - FULL SAMPLE (OLS, with Dx): forest plot on standardized beta, EAA only
#   Input: aggregate results table (group-level statistics; NO participant data).
#   Reads sheet "PRIMARY_complete" of EAA_primary_OLS_Dx.xlsx.
#   Outputs: BForest_Full_std.{png,svg}.
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
WORKDIR  <- getwd()
EAA_FILE <- file.path(WORKDIR, "EAA_primary_OLS_Dx.xlsx")   # sheet PRIMARY_complete (aggregate results)
while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
if (!dir.exists(WORKDIR)) dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
cat("Output folder:", getwd(), "\n")
FONT <- "Helvetica"
C_FDR <- "#D62728"   # red - FDR-significant
C_NOM <- "#E8887C"   # pink - nominal
nonsig_color   <- "#C2C2C2"
family_headers <- c("First/Second Generation Clocks", "Next Generation Clocks", "Organs/Systems Clocks")
fam_levels     <- c("First/Second Generation", "Next Generation", "Organs/Systems")
family_base    <- c("First/Second Generation" = "#9C27B0", "Next Generation" = "#F57B00", "Organs/Systems" = "#009E73")

shades <- function(base, n, lite=0.55, dark=0.40) {
  if (n == 1) return(base)
  brgb <- col2rgb(base)[,1] / 255
  mixw <- function(a) grDevices::rgb(brgb[1]+(1-brgb[1])*a, brgb[2]+(1-brgb[2])*a, brgb[3]+(1-brgb[3])*a)
  mixb <- function(a) grDevices::rgb(brgb[1]*(1-a), brgb[2]*(1-a), brgb[3]*(1-a))
  ts <- seq(lite, -dark, length.out = n)
  sapply(ts, function(t) if (t >= 0) mixw(t) else mixb(-t))
}

# ---- EAA rows: PRIMARY_complete = Measure | beta | 95% CI (std) | P | P-FDR | Cohen's d | Partial R2 ----
raw <- read_excel(EAA_FILE, sheet = "PRIMARY_complete")
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

dFull <- raw %>%
  filter(!is.na(Measure), !(Measure %in% family_headers)) %>%
  mutate(
    beta_std  = as.numeric(beta_std),
    P_num     = as.numeric(str_extract(as.character(P), "^[0-9.]+")),
    P_FDR_num = suppressWarnings(as.numeric(str_extract(as.character(P_FDR),"^[0-9.]+"))),
    CI_L_std  = as.numeric(str_extract(CI, "(?<=\\[)[-0-9.]+")),
    CI_U_std  = as.numeric(str_extract(CI, "(?<=,\\s)[-0-9.]+")),
    sig_p     = P_FDR_num,
    Label     = str_replace(Measure, "^EAA ", ""),
    Label     = ifelse(Label %in% c("IntrinsicCapacity","DunedinPACE"), paste0("* ", Label), Label)
  ) %>%
  mutate(.l = pmin(CI_L_std, CI_U_std), .u = pmax(CI_L_std, CI_U_std),
         CI_L_std = .l, CI_U_std = .u) %>% select(-.l, -.u)

dFull$Family <- factor(dFull$Family, levels = fam_levels)
dFull <- dFull %>% arrange(Family)
dFull$ClockColor <- NA_character_
for (fam in levels(dFull$Family)) {
  idx <- which(dFull$Family == fam); if (!length(idx)) next
  dFull$ClockColor[idx] <- shades(family_base[[fam]], length(idx))
}
sig_color <- function(pfdr, pnom) {
  if (!is.na(pfdr) && pfdr < 0.05) return(C_FDR)
  if (!is.na(pnom) && pnom < 0.05) return(C_NOM)
  return(nonsig_color)
}
dFull$FillColA  <- mapply(sig_color, dFull$sig_p, dFull$P_num)
dFull$FillSolid <- dFull$FillColA
df <- dFull %>% mutate(Label = factor(Label, levels = rev(Label)))

# Automatic range from the data (points and standardized CI bounds)
.betas <- suppressWarnings(as.numeric(dFull$beta_std))
.cil   <- suppressWarnings(as.numeric(dFull$CI_L_std))
.ciu   <- suppressWarnings(as.numeric(dFull$CI_U_std))
.dmin  <- min(c(.betas, .cil), na.rm = TRUE)
.dmax  <- max(c(.betas, .ciu), na.rm = TRUE)
x_lo <- min(-0.20, floor((.dmin - 0.01) / 0.05) * 0.05)
x_hi <- max( 0.25, ceiling((.dmax + 0.01) / 0.05) * 0.05)
span <- x_hi - x_lo
sq_w <- span * 0.035
sq_x <- x_lo - (span * 0.040)

fFull <- ggplot(df, aes(y = Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey50") +
  geom_errorbar(aes(xmin = CI_L_std, xmax = CI_U_std, colour = FillColA), orientation = "y", width = 0.28, linewidth = 0.7) +
  geom_point(aes(x = beta_std, colour = FillColA), size = 2.0) +
  geom_tile(aes(x = sq_x, fill = ClockColor), width = sq_w, height = 0.68, colour = "black", linewidth = 0.1) +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_x_continuous(expression(paste("Standardized Effect (", beta[std], ")")), breaks = seq(x_lo, x_hi, 0.1)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
  labs(y = NULL) +
  theme_classic(base_size = 7) +
  theme(text = element_text(family = FONT),
        axis.text.y  = element_text(size = 5.5, colour = "black", margin = margin(r = 4)),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        axis.text.x  = element_text(size = 4.5, colour = "black"),
        axis.title.x = element_text(size = 5),
        axis.line.x  = element_line(linewidth = 0.15),
        axis.ticks.x = element_line(linewidth = 0.15),
        plot.margin  = margin(2, 2, 2, 20))

mmin <- function(x) x / 25.4
while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
.p1 <- file.path(getwd(), "BForest_Full_std.png"); if (file.exists(.p1)) try(file.remove(.p1), silent=TRUE)
ggsave(.p1, fFull, width = mmin(55), height = mmin(62), dpi = 450, bg = "white")
while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
.p2 <- file.path(getwd(), "BForest_Full_std.svg"); if (file.exists(.p2)) try(file.remove(.p2), silent=TRUE)
ggsave(.p2, fFull, width = mmin(55), height = mmin(62), device = svglite::svglite, bg = "white")
while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
cat("Saved forest (EAA only, with Dx) in:", getwd(), "\n")
