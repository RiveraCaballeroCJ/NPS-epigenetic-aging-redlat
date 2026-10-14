# ============================================================================
# Figure 2 - FULL SAMPLE (OLS, with Dx): chord diagram + legend, EAA only
#   Input: aggregate results table (group-level statistics; NO participant data).
#   Reads sheet "PRIMARY_complete" of EAA_primary_OLS_Dx.xlsx.
#   Outputs: B1_legend_scale.{png,svg}, B2_chord.{png,svg}.
# ============================================================================
if (!require("readxl")) install.packages("readxl")
if (!require("circlize")) install.packages("circlize")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("stringr")) install.packages("stringr")
if (!require("svglite")) install.packages("svglite")
library(readxl); library(circlize); library(dplyr)
library(ggplot2); library(stringr); library(svglite); library(grid)

# ---- USER CONFIG: set this to your own output folder ----
# Figures (B1/B2) are written here, and the input xlsx is expected in this folder.
WORKDIR <- "set/your/output/folder/here"
while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
if (!dir.exists(WORKDIR)) { ok <- dir.create(WORKDIR, recursive = TRUE, showWarnings = TRUE)
if (!ok) stop("Could not create WORKDIR: ", WORKDIR) }
setwd(WORKDIR)
.tf <- file.path(WORKDIR, ".write_test.tmp")
if (inherits(try({writeLines("ok", .tf); file.remove(.tf)}, silent = TRUE), "try-error"))
  stop("No write permission in WORKDIR: ", WORKDIR)
cat("Output folder ready:", getwd(), "\n")

EAA_FILE <- file.path(WORKDIR, "EAA_primary_OLS_Dx.xlsx")   # sheet PRIMARY_complete
FONT <- "Helvetica"
nonsig_color   <- "#C2C2C2"; nonsig_ribbon  <- "#D3D3D3"
family_headers <- c("First/Second Generation Clocks","Next Generation Clocks","Organs/Systems Clocks")
fam_levels <- c("First/Second Generation", "Next Generation", "Organs/Systems")
fam_label  <- c("First/Second Generation"="First/Second Generation","Next Generation"="Next Generation",
                "Organs/Systems"="Organs/Systems")
family_base <- c("First/Second Generation"="#9C27B0","Next Generation"="#F57B00","Organs/Systems"="#009E73")

shades <- function(base, n, lite=0.55, dark=0.40) {
  if (n == 1) return(base)
  brgb <- col2rgb(base)[,1]/255
  mixw <- function(a) grDevices::rgb(brgb[1]+(1-brgb[1])*a, brgb[2]+(1-brgb[2])*a, brgb[3]+(1-brgb[3])*a)
  mixb <- function(a) grDevices::rgb(brgb[1]*(1-a), brgb[2]*(1-a), brgb[3]*(1-a))
  ts <- seq(lite, -dark, length.out = n)
  sapply(ts, function(t) if (t >= 0) mixw(t) else mixb(-t))
}

# EAA - PRIMARY_complete = Measure | beta | 95% CI (std) | P | P-FDR | Cohen's d | Partial R2
raw <- read_excel(EAA_FILE, sheet = "PRIMARY_complete")
raw <- raw[, 1:7]
names(raw) <- c("Measure","beta_std","CI","P","P_FDR","Cohen_d","R2")

cf <- NA; fv <- character(nrow(raw))
for (i in seq_len(nrow(raw))) {
  val <- raw$Measure[i]
  if (!is.na(val) && val %in% family_headers) {
    cf <- dplyr::case_when(
      val=="First/Second Generation Clocks" ~ "First/Second Generation",
      val=="Next Generation Clocks"         ~ "Next Generation",
      val=="Organs/Systems Clocks"          ~ "Organs/Systems")
  }
  fv[i] <- cf
}
raw$Family <- fv
dat <- raw %>%
  filter(!is.na(Measure), !(Measure %in% family_headers)) %>%
  mutate(beta_std=as.numeric(beta_std),
         P_num=as.numeric(str_extract(as.character(P),"^[0-9.]+")),
         P_FDR_num=suppressWarnings(as.numeric(str_extract(as.character(P_FDR),"^[0-9.]+"))),
         CI_L=as.numeric(str_extract(CI,"(?<=\\[)[-0-9.]+")),
         CI_U=as.numeric(str_extract(CI,"(?<=,\\s)[-0-9.]+"))) %>%
  mutate(sig_p=P_FDR_num,
         Label=str_replace(Measure, "^EAA ", ""))
dat$Family <- factor(dat$Family, levels = fam_levels)
dat <- dat %>% arrange(Family)
dat$ClockColor <- NA_character_
for (fam in fam_levels) {
  idx <- which(dat$Family == fam); if (!length(idx)) next
  dat$ClockColor[idx] <- shades(family_base[[fam]], length(idx))
}
dat$is_sig <- (!is.na(dat$sig_p) & dat$sig_p < 0.05) | (!is.na(dat$P_num) & dat$P_num < 0.05)
nonsig_ribbon <- "#E8E8E8"
dat$LinkColor <- ifelse(dat$is_sig, dat$ClockColor, nonsig_ribbon)

# ---- legend ----
fams_present <- fam_levels[sapply(fam_levels, function(f) any(dat$Family==f))]
NSEG<-24; bar_w<-1.15; bar_h<-0.21; row_gap<-0.62
seg_df<-list(); lab_df<-list(); frame_df<-list(); yrow<-0
for (fam in fams_present) {
  base <- family_base[[fam]]
  cols <- shades(base, NSEG)
  for (k in seq_len(NSEG)) { xa<-(k-1)/NSEG*bar_w; xb<-k/NSEG*bar_w
  seg_df[[length(seg_df)+1]]<-data.frame(xmin=xa,xmax=xb,ymin=yrow,ymax=yrow+bar_h,col=cols[k],stringsAsFactors=FALSE) }
  frame_df[[length(frame_df)+1]]<-data.frame(xmin=0,xmax=bar_w,ymin=yrow,ymax=yrow+bar_h)
  lab_df[[length(lab_df)+1]]<-data.frame(x=0,y=yrow-0.22,label=fam_label[fam],stringsAsFactors=FALSE)
  yrow<-yrow-row_gap
}
seg<-do.call(rbind,seg_df); labs<-do.call(rbind,lab_df); frames<-do.call(rbind,frame_df); ytop<-max(seg$ymax)
C_FDR<-"#D62728"; C_NOM<-"#E8887C"; C_NS<-"#C2C2C2"
ns_base<-min(labs$y)-1.15
sig3<-data.frame(col=c(C_FDR,C_NOM,C_NS),
                 label=c("italic(P)*'-FDR < 0.05'","italic(P)*'-value < 0.05'","italic(P)*' > 0.05'"),stringsAsFactors=FALSE)
sig3$sy<-ns_base-(seq_len(nrow(sig3))-1)*(bar_h+0.14)
nps_note_y<-min(labs$y)-0.55; nps_title_y<-ns_base+bar_h+0.20
clock_legend <- ggplot() +
  geom_rect(data=seg, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=col,colour=col), linewidth=0.3) +
  geom_rect(data=frames, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill=NA, colour="black", linewidth=0.2) +
  scale_fill_identity() + scale_colour_identity() +
  geom_text(data=labs, aes(x=x,y=y,label=label), hjust=0, vjust=1, size=6/.pt, family=FONT) +
  annotate("text", x=0, y=ytop+0.75, label="Aging clock families", hjust=0, size=7/.pt, family=FONT) +
  annotate("text", x=0, y=ytop+0.40, label="Epigenetic age acceleration (EAA)", hjust=0, size=5/.pt, family=FONT) +
  annotate("text", x=0, y=nps_note_y, label="* Native scale (index / pace ratio), not EAA", hjust=0, size=4.5/.pt, family=FONT) +
  annotate("text", x=0, y=nps_title_y, label="NPS effect", hjust=0, size=6.5/.pt, family=FONT) +
  geom_rect(data=sig3, aes(xmin=0, xmax=bar_h, ymin=sy, ymax=sy+bar_h, fill=col), colour="black", linewidth=0.3) +
  geom_text(data=sig3, aes(x=bar_h+0.12, y=sy+bar_h/2, label=label), parse=TRUE, hjust=0, size=5/.pt, family=FONT) +
  scale_x_continuous(limits=c(-0.1,3.6)) + scale_y_continuous(limits=c(min(sig3$sy)-0.3, ytop+1.1)) +
  theme_void() + theme(plot.margin=margin(2,4,2,2))

make_chord <- function(d, fam_text = FALSE) {
  circos.clear(); d <- d[!is.na(d$beta_std), ]
  mat <- matrix(pmax(abs(d$beta_std), 0.01), nrow = 1)
  rownames(mat) <- "NPS"; colnames(mat) <- d$Label
  all_sectors <- c("NPS", d$Label); gaps <- setNames(rep(1.5, length(all_sectors)), all_sectors); gaps["NPS"] <- 8.5
  circos.par(gap.after=gaps, start.degree=17.5, points.overflow.warning=FALSE,
             canvas.xlim=c(-1.35,1.35), canvas.ylim=c(-1.35,1.35), circle.margin=c(0.30,0.30,0.30,0.30))
  grid.col <- c("NPS"="#F0E442", setNames(d$ClockColor, d$Label))
  chordDiagram(mat, grid.col=grid.col, col=d$LinkColor, directional=0, transparency=0,
               annotationTrack="grid", link.border=NA, reduce=0, preAllocateTracks=list(track.height=0.10))
  circos.trackPlotRegion(track.index=1, bg.border=NA, panel.fun=function(x,y){
    sn<-get.cell.meta.data("sector.index"); xl<-get.cell.meta.data("xlim"); yl<-get.cell.meta.data("ylim")
    if (sn=="NPS") circos.text(mean(xl), yl[2]+0.05, "NPS", facing="inside", niceFacing=TRUE, cex=0.58, font=1)
    else { lbl <- if (sn %in% c("IntrinsicCapacity","DunedinPACE")) paste0("* ", sn) else sn
    circos.text(mean(xl), yl[2]+0.10, lbl, facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=0.45, font=1) }
  })
  if (fam_text) for (fam in levels(d$Family)) {
    secs <- as.character(d$Label[d$Family==fam]); if (!length(secs)) next
    highlight.sector(secs, track.index=1, col=NA, text=fam_label[[fam]], cex=0.50, text.col="black",
                     facing="bending.outside", niceFacing=TRUE) }
  circos.clear()
}

PH <- 100
draw_left <- function() { grid.newpage(); print(clock_legend, vp=viewport(width=1, height=1)) }
close_all_dev <- function() { while (!is.null(dev.list())) try(dev.off(), silent = TRUE) }
save_png <- function(fname, w_mm, h_mm, draw_fun) {
  close_all_dev(); path <- file.path(getwd(), fname); if (file.exists(path)) try(file.remove(path), silent=TRUE)
  grDevices::png(filename=path, width=w_mm/25.4, height=h_mm/25.4, units="in", res=450, pointsize=12)
  on.exit(try(dev.off(), silent=TRUE), add=TRUE); draw_fun() }
save_svg <- function(fname, w_mm, h_mm, draw_fun) {
  close_all_dev(); path <- file.path(getwd(), fname); if (file.exists(path)) try(file.remove(path), silent=TRUE)
  svglite::svglite(filename=path, width=w_mm/25.4, height=h_mm/25.4, pointsize=12)
  on.exit(try(dev.off(), silent=TRUE), add=TRUE); draw_fun() }
save_png("B1_legend_scale.png", 50, PH, draw_left)
save_svg("B1_legend_scale.svg", 50, PH, draw_left)
circos.clear(); save_png("B2_chord.png", 80, PH, function(){ circos.clear(); make_chord(dat) })
circos.clear(); save_svg("B2_chord.svg", 80, PH, function(){ circos.clear(); make_chord(dat) })
close_all_dev(); circos.clear()
cat("Saved B1 and B2 (EAA only, with Dx) in:", getwd(), "\n")
