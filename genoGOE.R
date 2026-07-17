# JMDplots/genoGEO.R
# Plots for the paper: Evolutionary oxidation of proteins in Earth's history
# 20240328 First JMDplots commit: eukaryotic gene age groups
# 20240409 Add Rubisco plots
# 20240528 Analyze methanogen genomes
# 20240803 Compare Rubisco proteins and unrelated genomes (stability_comparison)
# 20241224 Add Zc and stability diagram for S-cycling genomes
# 20250325 Add ancestral nitrogenase
# 20250625 Put all ancestral proteins in one plot and add thioredoxin and IPMDH
# 20250626 Put stability diagrams in one plot
# 20250627 Use average affinity instead of average rank of affinity for groupwise relative stability
# 20250903 Add Zc range diagram
# 20260313 Add BacDive analysis
# 20260319 Move to GitHub (jedick/genoGOE)

# Figure 1: Ranges of carbon oxidation state for organic compounds, amino acids, and proteins
genoGOE_1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_1.pdf", width = 11, height = 6)
  par(mar = c(0, 0, 0, 0))

  # Start with an empty plot
  plot(c(0, 10), c(0, 10), type = "n", axes = FALSE, xlab = "", ylab = "")
  # y-limits
  ybottom <- 0.5
  ytop <- 9.5
  # Get gradient colors
  col <- plotrix::smoothColors("2", 100, "4")
  # Loop over x-left values
  xs <- c(0, 3, 6)
  for(xleft in xs) {
    # Draw a smooth gradient across the bar
    plotrix::gradient.rect(xleft, ybottom, xleft + 0.5, ytop, col = col, border = NA, gradient = "y")
  }

  # Read protein data
  aa <- read.csv("UniProt/UP000000625_83333.csv.xz")
  Zc <- canprot::Zc(aa)
  # Print median Zc
  print(paste("Median Zc for E. coli proteins:", round(median(Zc), 3)))

  # Get 90% interval
  q90 <- quantile(Zc, probs = c(0.05, 0.95))
  # Get elemental formulas for reduced and oxidized proteins
  ilow <- which.min(abs(Zc - q90[1]))
  pf_low <- CHNOSZ::as.chemical.formula(CHNOSZ::protein.formula(aa[ilow, ]))
  ihigh <- which.min(abs(Zc - q90[2]))
  pf_high <- CHNOSZ::as.chemical.formula(CHNOSZ::protein.formula(aa[ihigh, ]))

  # Add labels

  text(xs[1], ybottom, "-4", cex = 2, font = 2, adj = c(0, 1.5))
  text(xs[1], ytop, "+4", cex = 2, font = 2, adj = c(0, -0.5))
  text(xs[1] + 0.6, ybottom, CHNOSZ::expr.species("CH4"), cex = 2, adj = 0)
  text(xs[1] + 0.6, ytop, CHNOSZ::expr.species("CO2"), cex = 2, adj = 0)

  text(xs[2], ybottom, "-1", cex = 2, font = 2, adj = c(0, 1.5))
  text(xs[2], ytop, "+1", cex = 2, font = 2, adj = c(0, -0.5))
  text(xs[2] + 0.6, ybottom, CHNOSZ::expr.species("C6H13NO2"), cex = 2, adj = 0)
  text(xs[2] + 0.6, ybottom + 0.8, "Leucine", cex = 2, adj = 0)
  text(xs[2] + 0.6, ytop, CHNOSZ::expr.species("C4H7NO4"), cex = 2, adj = 0)
  text(xs[2] + 0.6, ytop - 1, "Aspartic\nacid", cex = 2, adj = 0)

  text(xs[3], ybottom, round(q90[1], 3), cex = 2, font = 2, adj = c(0, 1.5))
  text(xs[3], ytop, round(q90[2], 3), cex = 2, font = 2, adj = c(0, -0.5))
  text(xs[3] + 1.1, ybottom, CHNOSZ::expr.species(pf_low), cex = 2, adj = 0)
  text(xs[3] + 1.1, ytop, CHNOSZ::expr.species(pf_high), cex = 2, adj = 0)

  # Add lines

  # Figure out y positions for amino acids range
  aa_range <- predict(lm(data.frame(y = c(ybottom, ytop), x = c(-4, 4))), data.frame(x = c(-1, 1)))
  lines(c(xs[1] + 0.53, xs[2] - 0.03), c(aa_range[1], ybottom + 0.02), lwd = 3, col = 2)
  lines(c(xs[1] + 0.53, xs[2] - 0.03), c(aa_range[2], ytop - 0.02), lwd = 3, col = 4)

  # Figure out y positions for protein range
  p_range <- predict(lm(data.frame(y = c(ybottom, ytop), x = c(-1, 1))), data.frame(x = c(q90[1], q90[2])))
  lines(c(xs[2] + 0.53, xs[3] - 0.03), c(p_range[1], ybottom + 0.02), lwd = 3, col = 2)
  lines(c(xs[2] + 0.53, xs[3] - 0.03), c(p_range[2], ytop - 0.02), lwd = 3, col = 4)

  # Add arrows and text for E. coli
  arrows(xs[3] + 0.7, ybottom + 0.02, xs[3] + 0.7, ytop - 0.02, code = 3, length = 0.2, angle = 35, lwd = 3)
  range_txt <- "90% of\nproteins\nin\nare in this\nrange"
  text(xs[3] + 1.1, mean(c(ybottom, ytop)), range_txt, cex = 2, adj = 0)
  E_coli_txt <- "    E. coli"
  text(xs[3] + 1.1, mean(c(ybottom, ytop)), E_coli_txt, cex = 2, adj = 0, font = 3)

  if(pdf) dev.off()
}

# Figure 2: Correlation of GC content and Zc with features from BacDive database
genoGOE_2 <- function(pdf = FALSE) {
  # Read data
  df <- read.csv("bacdive/cleaned_data_with_Zc.csv", check.names = FALSE)

  # Binary classification of Oxygen tolerance (Anaerobe or Non-anaerobe)
  is_anaerobe <- grepl("anaerobe", df$`oxygen_tolerance.Oxygen tolerance`, ignore.case = TRUE)
  df$`Oxygen tolerance` <- ifelse(is_anaerobe, "Anaerobe", "Non-anaerobe")

  # Keep only rows with all variables needed for plotting
  plot_df <- subset(
    df,
    !is.na(`culture_temp.Temperature`) &
      !is.na(`GC_content.GC-content`) &
      !is.na(Zc)
  )

  # Scatter: GC content vs temperature with linear fits by Oxygen tolerance
  p_gc <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = `culture_temp.Temperature`,
      y = `GC_content.GC-content`,
      color = `Oxygen tolerance`
    )
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Culture temperature (°C)",
      y = "GC content (%)",
      color = "Oxygen tolerance",
      title = "A"
    ) +
    ggplot2::scale_color_manual(values = c("Anaerobe" = "2", "Non-anaerobe" = "4")) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # Scatter: Zc vs temperature with linear fits by Oxygen tolerance
  p_zc <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = `culture_temp.Temperature`,
      y = Zc,
      color = `Oxygen tolerance`
    )
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Culture temperature (°C)",
      y = "Zc",
      color = "Oxygen tolerance",
      title = "B"
    ) +
    ggplot2::scale_color_manual(values = c("Anaerobe" = "2", "Non-anaerobe" = "4")) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # Arrange plots in a single row with shared legend
  combined_plot <- p_gc + p_zc + patchwork::plot_layout(nrow = 1, guides = "collect")

  combined_plot <- combined_plot & ggplot2::theme(legend.position = "bottom")

  if(pdf) file <- "Figure_2.pdf" else file <- "Figure_2.png"
  ggplot2::ggsave(file, combined_plot, width = 10, height = 4, dpi = 300)

}

# Figure 3: Genome-wide differences of oxidation state between two lineages of methanogens
genoGOE_3 <- function(pdf = FALSE, panel = NULL) {

  if(is.null(panel)) {
    if(pdf) pdf("Figure_3.pdf", width = 8, height = 6)
    mat <- matrix(c(1,2,3, 1,2,4, 5,5,5), nrow = 3, byrow = TRUE)
    layout(mat, heights = c(1, 1, 2))
  }
  panels <- if(is.null(panel)) LETTERS[1:4] else panel
  opar <- par(mgp = c(2.8, 1, 0), mar = c(5.1, 4.1, 2.1, 2.1))

  # Read methanogen genomes information
  mgfile <- "GTDB/methanogen_genomes.csv"
  mg <- read.csv(mgfile)
  # Genomes in Halobacteriota
  Halo <- mg$Genome[mg$Methanogen_class == "II"]
  # Genomes in Methanobacteriota
  Methano <- mg$Genome[mg$Methanogen_class == "I"]

  # Panels A-B: Zc and GC of marker genes 20240528

  # Read data for GTDB marker genes
  markerfile <- "GTDB/ar53_msa_marker_info_r220_XHZ+06.csv"
  markerdat <- read.csv(markerfile)
  markerid <- sapply(strsplit(markerdat$Marker.Id, "_"), "[", 2)
  methanogendir <- "methanogen"

  # NULL values for variables used in subset()
  protein <- gene <- NULL

  get_Zc <- function(phylum = "Halo") {
    files <- file.path(methanogendir, "marker", "faa", paste0(markerid, ".faa.xz"))
    aalist <- lapply(files, function(file) {
      aa <- suppressMessages(canprot::read_fasta(file))
      fivenum(canprot::Zc(subset(aa, protein %in% get(phylum))))
    })
    do.call(rbind, aalist)
  }

  get_GC <- function(phylum = "Halo") {
    files <- file.path(methanogendir, "marker", "fna", paste0(markerid, ".fna.xz"))
    nalist <- lapply(files, function(file) {
      na <- suppressMessages(canprot::read_fasta(file, molecule = "DNA"))
      na <- subset(na, gene %in% get(phylum))
      GC <- (na$G + na$C) / (na$G + na$C + na$A + na$T)
      fivenum(GC)
    })
    do.call(rbind, nalist)
  }

  # Get Zc for species in each phylum
  Zc_Halo <- get_Zc("Halo")
  Zc_Methano <- get_Zc("Methano")

  # Order by median Zc of Methanobacteriota
  iord <- order(Zc_Methano[, 3])
  Zc_Halo <- Zc_Halo[iord, ]
  Zc_Methano <- Zc_Methano[iord, ]

  if("A" %in% panels) {
    # Plot IQR of Zc
    plot(c(1, 53), c(-0.28, -0.04), xlab = "Marker gene", ylab = quote("Protein"~italic(Z)[C]), type = "n")
    for(i in 1:53) {
      lines(c(i, i) - 0.1, Zc_Methano[i, c(2, 4)], col = 2)
      lines(c(i, i) + 0.1, Zc_Halo[i, c(2, 4)], col = 4)
    }
    # Add legend for Class I and II methanogens
    legend("bottomright", "Class I", lty = 1, col = 2, bty = "n")
    legend("topleft", "Class II", lty = 1, col = 4, bty = "n")
    if(is.null(panel)) label.figure("A", font = 2, cex = 1.6)
    # Calculate p-value 20250304
    # Use median value in each group (3rd column) and paired observations
    p <- t.test(Zc_Halo[, 3], Zc_Methano[, 3], paired = TRUE)$p.value
    ptext <- bquote(italic(p) == .(signif(p, 2)))
    text(5, par("usr")[3], ptext, adj = c(0, -0.5))
  }

  # Get GC for species in each phylum
  GC_Halo <- get_GC("Halo")
  GC_Methano <- get_GC("Methano")
  GC_Halo <- GC_Halo[iord, ]
  GC_Methano <- GC_Methano[iord, ]

  if("B" %in% panels) {
    # Plot IQR of GC
    plot(c(1, 53), c(0.25, 0.65), xlab = "Marker gene", ylab = "GC content", type = "n")
    for(i in 1:53) {
      lines(c(i, i) - 0.1, GC_Methano[i, c(2, 4)], col = 2)
      lines(c(i, i) + 0.1, GC_Halo[i, c(2, 4)], col = 4)
    }
    if(is.null(panel)) label.figure("B", font = 2, cex = 1.6)
    # Calculate p-value 20250304
    # Use median value in each group (3rd column) and paired observations
    p <- t.test(GC_Halo[, 3], GC_Methano[, 3], paired = TRUE)$p.value
    ptext <- bquote(italic(p) == .(signif(p, 2)))
    text(5, par("usr")[3], ptext, adj = c(0, -0.5))
  }

  # Panel C: Delta Zc for marker genes

  if("C" %in% panels) {

    # If plotting only this panel, only make the abundance plot
    if(is.null(panel)) {
      # Plot Delta Zc vs Delta GC
      par(mar = c(4.1, 4.1, 1.1, 2.1))
      Delta_Zc <- na.omit(Zc_Halo[, 3] - Zc_Methano[, 3])
      Delta_GC <- na.omit(GC_Halo[, 3] - GC_Methano[, 3])
      plot(Delta_GC, Delta_Zc, xlab = quote(Delta*"GC"),
        ylab = quote(Delta*italic(Z)[C]~"(Class II - Class I)                                                         "),
        pch = 19, col = adjustcolor(1, alpha.f = 0.5), xpd = NA)
      # Calculate linear fit
      mylm <- lm(Delta_Zc ~ Delta_GC)
      x <- range(Delta_GC)
      y <- predict.lm(mylm, data.frame(Delta_GC = x))
      # Plot linear fit and show R2
      lines(x, y, lty = 2, lwd = 1.5, col = 8)
      R2 <- summary(mylm)$r.squared
      R2_txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 2, format = "f")))
      legend("topleft", legend = R2_txt, bty = "n", inset = c(-0.05, 0))
      if(is.null(panel)) label.figure("C", font = 2, cex = 1.6, yfrac = 0.9)
    }

    ylab <- if(is.null(panel)) "" else quote(Delta*italic(Z)[C]~"(Class II - Class I)")
    # Plot Delta Zc vs log10 protein abundance in M. maripaludis 20240531
    Delta_Zc <- Zc_Halo[, 3] - Zc_Methano[, 3]
    abundance <- markerdat$Redundant.Peptides / markerdat$MW
    log10a <- log10(abundance)
    plot(log10a, Delta_Zc, xlab = quote(log[10]~"protein abundance in"~italic("M. maripaludis")), ylab = ylab, pch = 19, col = adjustcolor(1, alpha.f = 0.5))
    # Calculate linear fit
    mylm <- lm(Delta_Zc ~ log10a)
    x <- range(log10a)
    y <- predict.lm(mylm, data.frame(log10a = x))
    # Plot linear fit and show R2
    lines(x, y, lty = 2, lwd = 1.5, col = 8)
    # Add horizontal line at Delta ZC = 0
    abline(h = 0, lty = 3)
    R2 <- summary(mylm)$r.squared
    R2_txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 2, format = "f")))
    legend("topleft", legend = R2_txt, bty = "n", inset = c(-0.05, 0))

  }

  par(opar)

  # Panel D: Zc controlled for various factors 20240529

  if("D" %in% panels) {
   
    # Get values of Zc, GC, and Cost
    genomes <- mg$Genome
    values <- lapply(genomes, function(genome) {
      aa <- read.csv(file.path(methanogendir, "aa", paste0(genome, "_aa.csv.xz")))
      data.frame(
        Zc = canprot::Zc(aa),
        GC = aa$abbrv,
        Cost = canprot::Cost(aa)
      )
    })
    names(values) <- genomes

    # NULL values for variables used in subset()
    GC <- Cost <- NULL
    # Get mean Zc for segment (phylum x condition)
    get_mean_Zc <- function(phylum = "Halo", condition = "all") {
      genomes <- get(phylum)
      myval <- values[genomes]
      if(condition == "low_GC") myval <- lapply(myval, subset, GC < 0.34)
      if(condition == "mid_GC") myval <- lapply(myval, subset, GC >= 0.34 & GC <= 0.36)
      if(condition == "high_GC") myval <- lapply(myval, subset, GC > 0.36)
      if(condition == "low_Cost") myval <- lapply(myval, subset, Cost < 23)
      if(condition == "mid_Cost") myval <- lapply(myval, subset, Cost >= 23 & Cost <= 24)
      if(condition == "high_Cost") myval <- lapply(myval, subset, Cost > 25)
      print(paste("median", median(sapply(myval, nrow)), "proteins for", phylum, condition))
      sapply(sapply(myval, "[", "Zc"), "mean")
    }

    Zc <- data.frame(
      Methano_all = get_mean_Zc("Methano", "all"),
      Halo_all = get_mean_Zc("Halo", "all"),
      Methano_low_GC = get_mean_Zc("Methano", "low_GC"),
      Halo_low_GC = get_mean_Zc("Halo", "low_GC"),
      Methano_mid_GC = get_mean_Zc("Methano", "mid_GC"),
      Halo_mid_GC = get_mean_Zc("Halo", "mid_GC"),
      Methano_high_GC = get_mean_Zc("Methano", "high_GC"),
      Halo_high_GC = get_mean_Zc("Halo", "high_GC"),
      Methano_low_Cost = get_mean_Zc("Methano", "low_Cost"),
      Halo_low_Cost = get_mean_Zc("Halo", "low_Cost"),
      Methano_mid_Cost = get_mean_Zc("Methano", "mid_Cost"),
      Halo_mid_Cost = get_mean_Zc("Halo", "mid_Cost"),
      Methano_high_Cost = get_mean_Zc("Methano", "high_Cost"),
      Halo_high_Cost = get_mean_Zc("Halo", "high_Cost")
    )

    # Don't plot overall line
    what <- c(FALSE, TRUE, TRUE, TRUE)
    # Start beanplot with all proteins
    par(mar = c(3.1, 4.1, 4.1, 2.1))
    bp <- beanplot::beanplot(Zc[, 1:2], side = "both", col = list(c(2, 7, 2, 2), c(4, 3, 4, 4)), xlim = c(0.5, 7.5), what = what, names = "")
    # Add means for species in each phylum
    abline(h = bp$stats[1], col = 2, lty = 2)
    abline(h = bp$stats[2], col = 4, lty = 2)
    # Add labels for methanogen classes
    text(0.6, bp$stats[1] + 0.008, "Class I")
    text(1.6, bp$stats[2] + 0.008, "Class II")
    # Add beans for GC and Cost
    beanplot::beanplot(Zc[, 3:14], side = "both", col = list(c(2, 7, 2, 2), c(4, 3, 4, 4)), xlim = c(0.5, 7.5), what = what, names = character(6), add = TRUE, at = 2:7)
    mtext(quote("Protein"~italic(Z)[C]), 2, line = 2.8, cex = par("cex"))

    # Add group names
    axis(1, at = 2:4, labels = c("GC < 0.34", "0.34 < GC < 0.36", "GC > 0.36"), gap.axis = 0)
    axis(1, at = 5:7, labels = c("Cost < 23", "23 < Cost < 25", "Cost > 25"), gap.axis = 0)
    axis(3, at = c(1, 3, 6), labels = c("Entire genomes", "Control for GC content", "Control for metabolic cost"), tick = FALSE, font = 2)

    if(is.null(panel)) label.figure("D", font = 2, cex = 1.6, xfrac = 0.018)

  }

  if(pdf & is.null(panel)) dev.off()

}

# Figure 4: Carbon oxidation state of proteins in eukaryotic gene age groups
genoGOE_4 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_4.pdf", width = 7, height = 6)
  layout(matrix(1:2), heights = c(1.2, 1.7))

  # Gene ages from Liebeskind et al. (2016)
  datadir <- "LMM16"
  modeAges <- read.csv(file.path(datadir, "modeAges_names.csv"))
  # Read list of reference proteomes
  refprot <- read.csv(file.path(datadir, "reference_proteomes.csv"))
  # Read summed amino acid compositions for proteins in each modeAge in each organism 20231218
  aa <- read.csv(file.path(datadir, "modeAges_aa.csv"))

  ylab <- "Carbon oxidation state of proteins"
  ylim <- c(-0.18, -0.06)
  yOE <- -0.07
  ytick.at <- seq(-0.18, -0.06, 0.04)
  ylabels.at <- c(-0.18, -0.06)

  # Loop over lineages
  coltext <- c("blue1", "green4")
  col <- adjustcolor(coltext, 0.6)
  lineages <- c("Mammalia", "Saccharomyceta")
  for(j in seq_along(lineages)) {
    # Adjust margins and colors
    if(j == 1) {
      par(mar = c(2.5, 4.4, 2, 3))
    }
    if(j == 2) {
      par(mar = c(7.3, 4.4, 2.5, 3))
    }
    # Start plot
    plot(c(1, 7), ylim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", yaxs = "i", type = "n", font.lab = 2)
    axis(2, ytick.at, labels = FALSE)
    axis(2, ylabels.at, tick = FALSE)
    # Get lineage-average Zc and uncertainty for each modeAge
    ilineage <- which(modeAges$X8 %in% lineages[j])
    lineage_vals <- lapply(ilineage, function(i) {
      OSCODE <- refprot$OSCODE[i]
      myaa <- aa[aa$organism == OSCODE, ]
      vals <- canprot::Zc(myaa)
      # Remove Euk+Bacteria (non-phylogenetic age category) 20231210
      vals <- vals[-3]
      vals
    })
    # Get Zc, mean and SD
    Zc_mat <- do.call(rbind, lineage_vals)
    Zc_mean <- colMeans(Zc_mat, na.rm = TRUE)
    Zc_sd <- apply(Zc_mat, 2, stats::sd, na.rm = TRUE)

    # Add single lineage line and error bars
    mode_ages <- 1:ncol(Zc_mat)
    lines(mode_ages, Zc_mean, col = col[j], lwd = 2)
    arrows(mode_ages, Zc_mean - Zc_sd, mode_ages, Zc_mean + Zc_sd,
      angle = 90, code = 3, length = 0.05, col = col[j]
    )
    if(j == 1) inset <- c(-0.02, 0.37) else inset <- c(-0.02, 0.22)
    legend("topleft", paste(lineages[j], "genomes", sep = "\n"), bty = "n", inset = inset)
    # Line for GOE
    lines(c(2.1, 2.5), c(yOE, yOE), lwd = 4, col = 2)
    if(j == 1) {
      # Top axis: GOE and NOE
      mtext("GOE", side = 3, at = 2.3, line = 0.5, font = 2)
      mtext("NOE", side = 3, at = 5, line = 0.5, font = 2)
      # Line for NOE
      lines(c(4.5, 6.5), c(yOE, yOE), lwd = 4, col = 2)
      # Divergence times (TimeTree 5)
      Mya <- c(4250, "3500-2600", 1598, 1275, 743, 563, 180)
      # Put Mya labels on bottom of first panel ...
      Mya_axis <- 1
    } else {
      # ... or top of second panel
      Mya_axis <- 3
      Mya <- c(4250, "3500-2600", 1598, 1275, 642, 528, 523)
      # Line for NOE
      lines(c(4.5, 5.5), c(yOE, yOE), lwd = 4, col = 2)
      # Bottom axis: tick marks and names for shared ancestry
      axis(1, at = 1:7, labels = FALSE)
      labels <- c("Cellular organisms", "Eukaryotes + Archaea", "Eukaryota", "Opisthokonta")
      text(x = 1:4, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 45, adj = 1, xpd = TRUE)
      # Bottom axis: names for Mammalia and Saccharomyceta lineages
      dx <- 0.2
      labels_slash <- c("/", "/", "/")
      text(x = (5:7) - dx, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels_slash, srt = 45, adj = 1, xpd = TRUE)
      labels_Mammalia <- c("Eumetazoa  ", "Vertebrata  ", "Mammalia  ")
      text(x = (5:7) - dx, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels_Mammalia, srt = 45, adj = 1, xpd = TRUE, col = coltext[1])
      labels_Saccharomyceta <- c("Dikarya", "Ascomycota", "Saccharomyceta")
      text(x = (5:7) + dx, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels_Saccharomyceta, srt = 45, adj = 1, xpd = TRUE, col = coltext[2])
    }
    # Add labels for divergence times
    at <- seq_along(Mya)
    axis(Mya_axis, at, Mya)
    # Add plot label 20240803
    label.plot(LETTERS[j], font = 2, cex = 1.2, xfrac = 0.04)
  }
  # Outer axis labels
  mtext(ylab, side = 2, line = 3, adj = -0.68, font = 2)
  mtext("Gene\nage\n(Ma)", side = 4, line = 1.5, font = 2, las = 1, at = -0.0215, adj = 0.5)

  if(pdf) dev.off()

}

# Figure 5: Carbon oxidation state of reconstructed ancestral sequences and extant proteins 20250625
genoGOE_5 <- function(pdf = FALSE) {

  # Function to add Zc labels with red/blue colors 20250625
  label_y_axis <- function() {
    # Add tick labels with red/blue colors 20250625
    axis(2, at = seq(-0.14, -0.12, 0.02), col.axis = 4)
    axis(2, at = -0.16)
    axis(2, at = seq(-0.28, -0.18, 0.02), col.axis = 2)
    # Add axis label
    mtext(quote(italic(Z)[C]), side = 2, line = 3.5, las = 0, cex = par("cex"))
  }

  # Plot Zc of rubisco from Kacar et al. (2017)  20240407
  plot_Rubisco <- function(ylim = c(-0.20, -0.12)) {
    # Read amino acid compositions
    fasta_file <- "KHAB17/rubisco.fasta"
    aa <- canprot::read_fasta(fasta_file)
    # Assign protein names
    aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

    # Get point locations
    xs <- 1:6
    ys <- canprot::Zc(aa)
    # Start plot
    xlab <- "Ancestral sequences (older to younger)"
    plot(xs, ys, type = "n", xaxt = "n", xlab = xlab, ylab = "", ylim = ylim, yaxt = "n")
    label_y_axis()
    # Plot main branch (excluding Anc I/III')
    lines(xs[-3], ys[-3], type = "b", pch = 19)
    # Add point for Anc I/III'
    points(xs[3], ys[3], pch = 19, col = 8)
    axis(1, at = 1:6, aa$protein)
    abline(v = 3.5, lty = 2, lwd = 2)
    text(2.9, -0.13, "GOE\n(estimated)")
  }


  # Plot Zc of IPMDH from Cui et al. (2025)  20250407
  plot_IPMDH <- function(ylim = c(-0.20, -0.12)) {
    aa <- canprot::read_fasta("CDY+25/IPMDH.fasta")
    Zc <- canprot::Zc(aa)
    # Ages from Table 1 of Cui et al., 2025
    ages <- c(
      2980, 2960, 2910, 2590,
      2360, 2160, 2140, 1570,
      1200, 932, 624, 0
    )
    # Uncertainties from Table 1 of Cui et al., 2025
    uncertainty <- c(
      200, 200, 200, 210,
      210, 220, 220, 250,
      270, 294, 318, 0
    )
    pch <- rep(19, length(Zc))
    pch[length(pch)] <- 1
    plot(ages / 1000, Zc, xlim = c(3.2, 0), xlab = "Age (Ga)", ylab = "", type = "b", ylim = ylim, yaxt = "n", pch = pch)
    # Add horizontal age uncertainty bars (age +/- uncertainty)
    i_err <- uncertainty > 0
    arrows(
      (ages[i_err] - uncertainty[i_err]) / 1000, Zc[i_err],
      (ages[i_err] + uncertainty[i_err]) / 1000, Zc[i_err],
      angle = 90, code = 3, length = 0.04
    )
    label_y_axis()
    # Lines for GOE and NOE
    yOE <- -0.12
    lines(c(2.5, 2.2), c(yOE, yOE), lwd = 4, col = 2)
    text(2.35, yOE - 0.005, "GOE")
    lines(c(0.8, 0.54), c(yOE, yOE), lwd = 4, col = 2)
    text(0.67, yOE - 0.005, "NOE")
  }

  # Plot Zc for ancestral thioredoxins from Perez-Jimenez et al. (2011)  20250625
  plot_thioredoxin <- function(ylim = c(-0.28, -0.18)) {
    # Read data file with ages and PDB IDs from Del Galdo et al. (2019)
    dat <- read.csv("PIZ+11/DAAD19.csv")
    # Read amino acid compositions
    aa <- canprot::read_fasta("PIZ+11/thioredoxin.fasta")
    # Calculate Zc
    Zc <- canprot::Zc(aa)
    # Setup plot
    xlim <- c(4.5, 0)
    plot(xlim, range(Zc), xlim = xlim, xlab = "Age (Ga)", ylab = "", type = "n", ylim = ylim, yaxt = "n")
    label_y_axis()
    # Helper to draw horizontal age uncertainty bars (skip time-zero points)
    add_age_error_bars <- function(i, yvals) {
      i_err <- i & !(dat$Age == 0 & dat$Min == 0 & dat$Max == 0)
      arrows(dat$Min[i_err], yvals[i_err], dat$Max[i_err], yvals[i_err],
        angle = 90, code = 3, length = 0.04
      )
    }
    # Add separate lines for each lineage
    iBac <- dat$Lineage == "Bacteria"
    add_age_error_bars(iBac, Zc)
    pch <- rep(19, sum(iBac))
    pch[length(pch)] <- 1
    lines(dat$Age[iBac], Zc[iBac], type = "b", pch = pch)
    text(3.5, -0.22, "Bacteria")
    iArcEuk <- dat$Lineage == "Arc-Euk"
    add_age_error_bars(iArcEuk, Zc)
    pch <- rep(15, sum(iArcEuk))
    pch[length(pch)] <- 0
    lines(dat$Age[iArcEuk], Zc[iArcEuk], type = "b", pch = pch)
    text(2.5, -0.252, "Archaea+Eukaryota")
    # Lines for GOE and NOE
    yOE <- -0.18
    lines(c(2.5, 2.2), c(yOE, yOE), lwd = 4, col = 2)
    text(2.35, yOE - 0.006, "GOE")
    lines(c(0.8, 0.54), c(yOE, yOE), lwd = 4, col = 2)
    text(0.67, yOE - 0.006, "NOE")
  }

  if(pdf) pdf("Figure_5.pdf", width = 9, height = 6)
  par(mfrow = c(2, 2))
  par(las = 1)
  par(mar = c(4.0, 5.0, 2.5, 1.0), mgp = c(2.5, 1, 0))
  plot_thioredoxin()
  title("Thioredoxin")
  label.figure("A", font = 2, cex = 1.5)
  plot_IPMDH()
  title("IPMDH")
  label.figure("B", font = 2, cex = 1.5)
  plot_Rubisco()
  title("Rubisco")
  label.figure("C", font = 2, cex = 1.5)
  plot_nitrogenase()
  title("Nitrogenase")
  label.figure("D", font = 2, cex = 1.5)
  if(pdf) dev.off()

}

# New version of Figure 5 with stacked Zc plots above O2 profile  20260715 jmd
genoGOE_5 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure_5.pdf", width = 8, height = 8)
  # Setup figure region
  par(mfrow = c(6, 1))
  par(mar = c(0, 4.1, 0, 2.1))
  # Make plots
  plot_phylostrata()
  plot_rubisco()
  plot_nitrogenase()
  plot_thioredoxin()
  plot_IPMDH()
  plot_oxygen()
  if(pdf) dev.off()
}

#########################################
### Supporting functions for Figure 5 ###
#########################################

plot_phylostrata <- function() {
  # Gene ages from Liebeskind et al. (2016)
  datadir <- "LMM16"
  modeAges <- read.csv(file.path(datadir, "modeAges_names.csv"))
  # Read list of reference proteomes
  refprot <- read.csv(file.path(datadir, "reference_proteomes.csv"))
  # Read summed amino acid compositions for proteins in each modeAge in each organism 20231218
  aa <- read.csv(file.path(datadir, "modeAges_aa.csv"))

  # Start plot
  ylim <- c(-0.20, -0.05)
  xlim <- c(4.5, 0)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  axis(1, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-0.20, -0.05, 0.05), labels = FALSE)
  axis(2, c(-0.20, -0.05), tick = FALSE, las = 1)

  # Lineage names, node ages, and colors
  lineages <- c("Mammalia", "Saccharomyceta")
  ages_Mya <- list(
    c(4250, mean(c(3500, 2600)), 1598, 1275, 743, 563, 180),
    c(4250, mean(c(3500, 2600)), 1598, 1275, 642, 528, 523)
  )
  coltext <- c("blue1", "green4")
  col <- adjustcolor(coltext, 0.6)

  # Loop over lineages
  for(j in seq_along(lineages)) {
    # Get lineage-average Zc and uncertainty for each modeAge
    ilineage <- which(modeAges$X8 %in% lineages[j])
    lineage_vals <- lapply(ilineage, function(i) {
      OSCODE <- refprot$OSCODE[i]
      myaa <- aa[aa$organism == OSCODE, ]
      vals <- canprot::Zc(myaa)
      # Remove Euk+Bacteria (non-phylogenetic age category) 20231210
      vals <- vals[-3]
      vals
    })
    # Get Zc, mean and SD
    Zc_mat <- do.call(rbind, lineage_vals)
    Zc_mean <- colMeans(Zc_mat, na.rm = TRUE)
    Zc_sd <- apply(Zc_mat, 2, stats::sd, na.rm = TRUE)

    # Add single lineage line and error bars
    ages_Gya <- ages_Mya[[j]] / 1000
    lines(ages_Gya, Zc_mean, col = col[j], type = "b", pch = substr(lineages[j], 1, 1))
    arrows(ages_Gya, Zc_mean - Zc_sd, ages_Gya, Zc_mean + Zc_sd,
      angle = 90, code = 3, length = 0.05, col = col[j]
    )
  }

  # Add legend
  legend("topright", lineages, pch = substr(lineages, 1, 1), col = coltext, bty = "n")
  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("Age (Ga)", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("Gene age groups", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("(Phylostrata: Liebeskind et al., 2016)", side = 3, line = -1.5, adj = 0.25, cex = par("cex") * 1.2)

}

# Plot Zc of rubisco from Kacar et al. (2017)  20240407
plot_rubisco <- function() {
  # Read amino acid compositions
  fasta_file <- "KHAB17/rubisco.fasta"
  aa <- canprot::read_fasta(fasta_file)
  # Assign protein names
  aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

  # Branch lengths from Fig. 4 of Kacar et al. (2017) (root node at zero)
  # Root node, Anc. I/II/III, Anc. I/III, Anc. I/III', Anc. I, Anc. IB, Anc. IAB, Nostoc-Anabaena split
  aa_sub <- c(0, 0.34, 0.59, 0.69, 1.63, 1.91, 1.95, 2.06)
  # Transform branch length to time (linear model)
  # Age of Rubisco ca. 3.8 Gya - Taylor-Kearney et al., 2024
  # Nostoc-Anabaena split median 1373 Mya, adjusted 1670 Mya - TimeTree accessed on 2026-07-15
  sub <- c(head(aa_sub, 1), tail(aa_sub, 1))
  age <- c(3.8, 1.373)
  sub2age <- lm(age ~ sub)
  # Predict ages from branch lengths for ancestral sequences
  ages <- predict(sub2age, newdata = data.frame(sub = head(tail(aa_sub, -1), -1)))
  # Make inverse model to get tick labels
  age2sub <- lm(sub ~ age)

  # Get Zc values
  Zc_vals <- canprot::Zc(aa)

  # Start plot
  ylim <- c(-0.20, -0.10)
  xlim <- c(4.5, 0)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  # Put labels at the Rubisco age and whole Gya after that
  at <- c(3.8, 3, 2, 1)
  labels <- round(predict(age2sub, newdata = data.frame(age = at)), 1)
  axis(1, at = at, labels = labels, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-0.20, -0.10, 0.05), labels = FALSE)
  axis(2, c(-0.20, -0.10), tick = FALSE, las = 1)

  # Plot main branch (excluding Anc I/III')
  lines(ages[-3], Zc_vals[-3], type = "b", pch = 19)
  # Add point for Anc I/III'
  points(ages[3], Zc_vals[3], pch = 19, col = 8)
  # Add text labels
  text(ages[1:4], Zc_vals[1:4] + 0.007, aa$protein[1:4])
  text(ages[5], Zc_vals[5] + 0.007, aa$protein[5], adj = 1)
  text(ages[6], Zc_vals[6] - 0.007, aa$protein[6], adj = 0)

  # Add legend
  legend("topright", legend = "Linear age model\nRelative ages only", bty = "n", text.font = 3, cex = 1.2)
  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("aa subs/site", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("Rubisco large subunit", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("(RAS data: Kacar et al., 2017)", side = 3, line = -1.5, adj = 0.30, cex = par("cex") * 1.2)
}

# Plot Zc for ancestral and modern nitrogenases from Cuevas Zuviría et al. (2025)  20260715
plot_nitrogenase <- function() {
  # Directory with FASTA files
  seq_dir <- "CDA+25"
  # Start of file name
  file_start <- "AGNifAlign105.ext-anc.alt"
  # Subunit names
  subunits <- c("D", "K", "H")
  # Node names (oldest to youngest)
  nodes <- c("1206_map", "1207_map", "1209_map", "1224_map", "Nif_Azotobacter_vinelandii")

  # Branch lengths from Fig. 3 of Cuevas Zuviría et al. (2025) (root node at zero)
  # Node substitutions (aa/site)
  # Root node, Anc1206, Anc1207, Anc1209, Anc1224, A. vinelandii
  aa_sub <- c(0, 0.23, 0.71, 0.93, 1.04, 1.76)

  # Transform branch length to time (linear model)
  # ca. 3.2 Gya for nitrogenase and 2.5 Gya for Anc1 - Cuevas Zuviría et al.
  sub <- c(0, 1.04)
  age <- c(3.2, 2.5)
  sub2age <- lm(age ~ sub)
  # Predict ages from branch lengths for ancestral and modern sequences
  ages <- predict(sub2age, newdata = data.frame(sub = tail(aa_sub, -1)))
  # Make inverse model to get tick labels
  age2sub <- lm(sub ~ age)

  # Start plot
  ylim <- c(-0.20, -0.12)
  xlim <- c(4.5, 0)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  # Put labels at the nitrogenase age and whole Gya after that
  at <- c(3.2, 3, 2)
  labels <- round(predict(age2sub, newdata = data.frame(age = at)), 1)
  axis(1, at = at, labels = labels, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-0.20, -0.12, 0.04), labels = FALSE)
  axis(2, c(-0.20, -0.12), tick = FALSE, las = 1)

  # Get Zc and nC at nodes for each subunit
  Zc_list <- nC_list <- vector("list", 3)
  for(i in 1:length(subunits)) {
    # Load amino acid composition from FASTA file
    fasta_file <- file.path(seq_dir, paste(file_start, subunits[i], "fasta", sep = "."))
    aa <- read_fasta(fasta_file)
    # Find the node names in the data frame
    iaa <- match(nodes, aa$protein)
    node_aa <- aa[iaa, ]
    # Calculate Zc and nC
    Zc_list[[i]] <- Zc(node_aa)
    nC_list[[i]] <- nC(node_aa)
    # Plot data points
    lines(ages, Zc_list[[i]], type = "b", pch = subunits[i])
  }

  ## Calculate mean Zc for DDKK (stoichiometry of Nif-I complex)
  #Zc_DDKK <- (Zc_list[[1]] * nC_list[[1]] + Zc_list[[2]] * nC_list[[2]]) / (nC_list[[1]] + nC_list[[2]])
  ## Plot points for DDKK
  #lines(ages, Zc_DDKK, type = "b", pch = 19)

  # Add legend
  legend("topright", legend = "Linear age model\nRelative ages only", bty = "n", text.font = 3, cex = 1.2)
  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("aa subs/site", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("Nitrogenase subunits", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("(RAS data: Cuevas Zuviría et al., 2025)", side = 3, line = -1.5, adj = 0.32, cex = par("cex") * 1.2)
}

# Plot Zc for ancestral thioredoxins from Perez-Jimenez et al. (2011)  20250625
plot_thioredoxin <- function() {
  # Read data file with ages and PDB IDs from Del Galdo et al. (2019)
  dat <- read.csv("PIZ+11/DAAD19.csv")
  # Read amino acid compositions
  aa <- canprot::read_fasta("PIZ+11/thioredoxin.fasta")
  # Calculate Zc
  Zc <- canprot::Zc(aa)

  # Start plot
  ylim <- c(-0.30, -0.15)
  xlim <- c(4.5, 0)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  axis(1, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-0.30, -0.15, 0.05), labels = FALSE)
  axis(2, c(-0.30, -0.15), tick = FALSE, las = 1)

  # Helper to draw horizontal age uncertainty bars (skip time-zero points)
  add_age_error_bars <- function(i, yvals) {
    i_err <- i & !(dat$Age == 0 & dat$Min == 0 & dat$Max == 0)
    arrows(dat$Min[i_err], yvals[i_err], dat$Max[i_err], yvals[i_err],
      angle = 90, code = 3, length = 0.04
    )
  }
  # Add separate lines for each lineage
  iBac <- dat$Lineage == "Bacteria"
  add_age_error_bars(iBac, Zc)
  pch <- rep(19, sum(iBac))
  pch[length(pch)] <- 1
  lines(dat$Age[iBac], Zc[iBac], type = "b", pch = pch)
  text(3.5, -0.22, "Bacteria")
  iArcEuk <- dat$Lineage == "Arc-Euk"
  add_age_error_bars(iArcEuk, Zc)
  pch <- rep(15, sum(iArcEuk))
  pch[length(pch)] <- 0
  lines(dat$Age[iArcEuk], Zc[iArcEuk], type = "b", pch = pch)
  text(2.5, -0.26, "Archaea+Eukaryota")

  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("Age (Ga)", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("Thioredoxin", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext(CHNOSZ::hyphen.in.pdf("(RAS data: Perez-Jimenez et al., 2011)"), side = 3, line = -1.5, adj = 0.19, cex = par("cex") * 1.2)
}

# Plot Zc of IPMDH from Cui et al. (2025)  20250407
plot_IPMDH <- function() {
  aa <- canprot::read_fasta("CDY+25/IPMDH.fasta")
  Zc <- canprot::Zc(aa)
  # Ages from Table 1 of Cui et al., 2025
  ages <- c(
    2980, 2960, 2910, 2590,
    2360, 2160, 2140, 1570,
    1200, 932, 624, 0
  )
  # Uncertainties from Table 1 of Cui et al., 2025
  uncertainty <- c(
    200, 200, 200, 210,
    210, 220, 220, 250,
    270, 294, 318, 0
  )
  pch <- rep(19, length(Zc))
  pch[length(pch)] <- 1

  # Start plot
  ylim <- c(-0.20, -0.12)
  xlim <- c(4.5, 0)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  axis(1, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-0.20, -0.12, 0.04), labels = FALSE)
  axis(2, c(-0.20, -0.12), tick = FALSE, las = 1)

  # Add points
  points(ages / 1000, Zc, type = "b", pch = pch)

  # Add horizontal age uncertainty bars (age +/- uncertainty)
  i_err <- uncertainty > 0
  arrows(
    (ages[i_err] - uncertainty[i_err]) / 1000, Zc[i_err],
    (ages[i_err] + uncertainty[i_err]) / 1000, Zc[i_err],
    angle = 90, code = 3, length = 0.04
  )

  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("Age (Ga)", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("IPMDH", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("(RAS: Cui et al., 2025)", side = 3, line = -1.5, adj = 0.10, cex = par("cex") * 1.2)
}

# Plot O2 curve from Lyons et al. (2024)  20260716
plot_oxygen <- function() {
  # Ribbon color
  col <- "#3E7CB1"

  # Start plot
  xlim <- c(4.5, 0)
  ylim <- c(-6, 0.2)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  axis(1, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-4, 0, 1), labels = FALSE)
  axis(2, c(-4, -2, 0), tick = FALSE, las = 1)
  # Ticks below breaks in scale
  axis(2, c(-5, -5.75), labels = c(-7, -13), las = 1)
  # Scale break marks
  text(par("usr")[1], -4.6, "~", xpd = NA, cex = 1.5, srt = 15)
  text(par("usr")[1], -5.5, "~", xpd = NA, cex = 1.5, srt = 15)

  # Archean
  lines(c(4, 3), c(-5.6, -5.6), col = col, lwd = 15)
  lines(c(3, 2.45), c(-5, -5), col = col, lwd = 15)

  # Plot smooth line with spline function
  smoothline <- function(x, y, x_vals, lwd = 10, lty = 1, rev = FALSE) {
    # x, y - x and y data
    # x_vals - x values for plotting
    spline_fun <- splinefun(x, y, method = "monoH.FC")
    y_vals <- spline_fun(x_vals)
    lines(x_vals, y_vals, col = col, lwd = lwd, lty = lty)
  }

  # GOE
  GOE_x <- c(2.431, 2.428, 2.416, 2.393, 2.373, 2.342, 2.312, 2.274, 2.236, 2.186, 2.15, 2.08, 2.021)
  GOE_y <- c(-5.114, -4.734, -4.152, -3.62, -3.316, -3.038, -2.81, -2.633, -2.506, -2.405, -2.354, -2.354, -2.354)
  smoothline(GOE_x, GOE_y, seq(2.43, 2.02, length.out = 20))

  # Overshoot
  OS_x <- c(2.177, 2.171, 2.163, 2.158, 2.15, 2.135, 2.107, 2.085, 2.075, 2.067, 2.054, 2.047, 2.034)
  OS_y <- c(-2, -1.62, -1.241, -0.861, -0.481, -0.101, 0, -0.354, -0.734, -1.139, -1.519, -1.899, -2.076)
  points(OS_x, OS_y, col = col, pch = 19)

  # Boring Billion
  #lines(c(2.02, 0.8), c(-2.35, -2.35), col = col, lwd = 25)
  # Shorten this thick line a bit for better alignment with adjoining (thinner) lines
  lines(c(1.97, 0.85), c(-2.35, -2.35), col = col, lwd = 25)

  # NOE
  NOE_x <- c(0.798, 0.793, 0.786, 0.766, 0.743, 0.733, 0.71, 0.687)
  NOE_y <- c(-2.127, -1.392, -1.038, -1.013, -1.19, -1.367, -1.544, -1.62)
  smoothline(NOE_x, NOE_y, seq(0.8, 0.69, length.out = 20))

  # Post-NOE
  lines(c(0.68, 0.42), c(-1.62, -1.62), col = col, lwd = 10)

  # POE
  POE_x <- c(0.422, 0.402, 0.387, 0.369, 0.321, 0.28, 0.232, 0.21, 0.144, 0.091, 0.033, 0.003)
  POE_y <- c(-1.595, -1.392, -1.089, -0.506, -0.127, 0.152, -0.051, -0.051, -0.152, 0.025, -0.076, 0.025)
  smoothline(POE_x, POE_y, seq(0.42, 0, length.out = 40))

  # Add axis labels and title
  mtext(quote(bolditalic(p)*bold(O[2])), side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2, padj = -2.3)
  mtext(quote(bold("("*log[10])), side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2, adj = 0.9, padj = 0)
  mtext("PAL)", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2, padj = 1.3)
  mtext("Age (Ga)", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("Oxygen", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("(Model: Lyons et al., 2024)", side = 3, line = -1.5, adj = 0.12, cex = par("cex") * 1.2)
}

#####################################
### End of functions for Figure 5 ###
#####################################

# Figure 6: Relative stabilities of proteins as a function of environmental variables
genoGOE_6 <- function(pdf = FALSE, panel = NULL) {

  if(is.null(panel)) {
    if(pdf) pdf("Figure_6.pdf", width = 12, height = 8)
    layout(matrix(c(1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3, 4,4,4,
                    rep(5, 9), rep(6, 9)), nrow = 2, byrow = TRUE))
    par(cex = 1)
  }
  panels <- if(is.null(panel)) LETTERS[1:7] else panel

  # Read amino acid compositions
  fasta_file <- "KHAB17/rubisco.fasta"
  aa <- canprot::read_fasta(fasta_file)
  # Assign protein names
  aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

  # Add proteins to CHNOSZ
  ip <- add.protein(aa, as.residue = TRUE)
  # Setup basis species and swap O2 for e- to make Eh-pH diagram
  basis("QEC+")
  swap.basis("O2", "e-")
  # Set plot resolution
  res <- 300
  
  # Panel A: Pairwise stability boundaries for Rubisco

  if("A" %in% panels) {

    # Loop over individual pairs
    for(pre in 1:3) {
      for(post in 4:6) {
        add <- TRUE
        if(pre == 1 & post == 4) add <- FALSE
        a <- affinity(pH = c(0, 14, res), Eh = c(-0.5, 0.8, res), iprotein = ip[c(pre, post)])
        d <- diagram(a, names = "", lty = 2, col = "#000000b0", add = add, limit.water = !add, fill.NA = "gray80", xlab = "pH", ylab = axis.label("Eh"))
        # Only label lines for reaction with Anc. I
        if(post == 4) {
          # Sort x values and get x and y values of boundary line
          order <- order(d$linesout[[1]])
          xs <- d$linesout[[1]][order]
          ys <- d$linesout[[2]][order]
          # Get a single value along the length of the line
          ilab <- floor(length(xs) * pre * 3 / 10)
          x <- xs[ilab]
          y <- ys[ilab]
          text(x, y - 0.03, aa$protein[pre], cex = 0.6)
          text(x, y + 0.03, aa$protein[post], cex = 0.6)
        }
      }
    }

    text(5.5, 0.67, CHNOSZ::hyphen.in.pdf("Higher affinity\nfor post-GOE protein\nin each pair"), cex = 0.8)
    text(4.8, -0.15, CHNOSZ::hyphen.in.pdf("Higher affinity for\npre-GOE protein in each pair"), cex = 0.8, srt = -37)
    title("Pairwise Rubiscos", font.main = 1)
    if(is.null(panel)) label.figure("A", cex = 1.5, font = 2, yfrac = 0.936)

  }

  if("B" %in% panels) {

    # Panel B: Groupwise stability boundaries for Rubisco

    # Calculate affinity of composition reactions for all proteins
    aout <- affinity(pH = c(0, 14, res), Eh = c(-0.5, 0.8, res), iprotein = ip)
    # Set up groups for affinity ranking:
    # 3 pre-GOE and 3 post-GOE proteins
    groups <- list(pre = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), post = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
    amean <- agg.affinity(aout, groups = groups)
    diagram(amean, lwd = 2, col = 4, names = "", xlab = "pH", ylab = axis.label("Eh"), balance = 1)
    text(6, -0.17, CHNOSZ::hyphen.in.pdf("Higher mean affinity\nfor pre-GOE proteins"), col = 4, font = 2, cex = 0.8, srt = -33)
    text(6.5, 0.1, CHNOSZ::hyphen.in.pdf("Higher mean affinity\nfor post-GOE proteins"), col = 4, font = 2, cex = 0.8, srt = -33)
    title("Groupwise Rubiscos", font.main = 1)
    if(is.null(panel)) label.figure("B", cex = 1.5, font = 2, yfrac = 0.936)

  }

  if("C" %in% panels) {

    # Panel D: Comparison between Rubiscos and methanogen and Nitrososphaeria genomes
    stability_comparison(res = res)

    # Label lines
    text(6.8, -61.4, "Rubiscos", cex = 0.8, adj = 0)
    text(6.5, -61.8, CHNOSZ::hyphen.in.pdf("Pre-GOE"), cex = 0.75, adj = 1, srt = 30)
    text(6.5, -60.7, CHNOSZ::hyphen.in.pdf("Post-GOE"), cex = 0.75, adj = 1, srt = 30)

    text(6.8, -67, "Methanogens", cex = 0.8, adj = 0)
    text(6, -67.6, "Class I", cex = 0.75)
    text(6, -66.7, "Class II", cex = 0.75)

    text(6.8, -69.4, "Nitrososphaeria", cex = 0.8, adj = 0, font = 3)
    text(6, -70.2, "Basal", cex = 0.75)
    text(6, -69.4, "Terrestrial", cex = 0.75)

    if(is.null(panel)) {
      title("   Rubiscos and all proteins in genomes", font.main = 1, xpd = NA)
      label.figure("C", cex = 1.5, font = 2, yfrac = 0.936)
    }

  }

  if("D" %in% panels) {

    # Add more arrows 20240812
    plot.new()
    opar <- par(xpd = NA)
    arrows(-0.3, 0.08, -0.3, 0.28, length = 0.2, lwd = 2, col = 2)
    arrows(-0.1, 0.22, -0.1, 0.42, length = 0.2, lwd = 2, col = 7)
    arrows(-0.3, 0.65, -0.3, 0.85, length = 0.2, lwd = 2, col = 4)
    text(0.05, 0.25, "Oxidation in\nmany lineages\naround GOE", adj =0)
    text(-0.15, 0.75, CHNOSZ::hyphen.in.pdf("Rubisco transitions\nat more oxidizing\nconditions"), adj = 0)
    par(opar)

  }

  # Evolutionary oxidation and relative stabilities for genomes with S-cycling genes 20241211
  # genoGOE/sulfur_genomes.R
  # 20241211 add Eh-pH affinity ranking
  # 20241223 convert to logfO2-pH

  # List genomes with single sulfur-cycling genes
  genomes <- list(
    dsrAB = c("GCA_002782605.1", "GCF_000517565.1", "GCA_002878135.1"),
    soxC = c("GCF_000153205.1", "GCF_000024725.1", "GCA_001914955.1", "GCA_002731275.1", 
      "GCA_002007425.1", "GCA_001780165.1", "GCA_002721445.1", "GCF_900129635.1", 
      "GCA_002162915.1", "GCA_002712885.1", "GCF_000484535.1", "GCF_000969705.1", 
      "GCF_002148795.1", "GCF_002514725.1", "GCF_900106035.1", "GCA_002687025.1", 
      "GCA_003222815.1", "GCF_900187885.1", "GCA_002712165.1", "GCA_002705185.1", 
      "GCA_003228115.1"),
    soxABXYZ = c("GCF_000021565.1", "GCF_900142435.1", "GCA_000830255.1", "GCF_000227215.1"),
    aprAB = c("GCA_001800245.1", "GCA_002898195.1", "GCA_002717185.1", "GCF_000328625.1", 
      "GCF_002252565.1", "GCA_001805205.1", "GCA_001784555.1", "GCA_001443375.1"),
    dmsA = c("GCF_000384115.1", "GCA_001593855.1", "GCF_000020005.1", "GCA_003242675.1", 
      "GCA_001771285.1", "GCA_002717245.1", "GCA_003223635.1", "GCF_001051235.1", 
      "GCA_001304035.1", "GCF_000772535.1", "GCA_002898895.1", "GCF_001049895.1", 
      "GCF_001860525.1", "GCA_001775395.1", "GCA_001775995.1", "GCA_001830835.1", 
      "GCF_000487995.1", "GCA_002747435.1", "GCA_002839495.1", "GCA_001515205.2", 
      "GCA_001742785.1", "GCA_001775755.1", "GCA_001768675.1"),
    mddA = c("GCA_003223145.1", "GCA_001872725.1", "GCF_000018105.1", "GCA_001447805.1", 
      "GCA_002400775.1", "GCA_001563325.1", "GCA_002746235.1", "GCA_001780825.1", 
      "GCF_000970205.1", "GCA_002746185.1", "GCA_002790835.1", "GCF_001886815.1", 
      "GCF_000192575.1", "GCA_002256595.1", "GCF_900111015.1", "GCA_001664505.1", 
      "GCA_002841995.1", "GCA_002699105.1", "GCF_002563855.1", "GCA_003219195.1"),
    dmdA = c("GCA_002707655.1", "GCF_900102465.1", "GCA_001800745.1", "GCF_001029505.1", 
      "GCA_002717565.1", "GCA_002701885.1", "GCA_002722565.1")
  )

  # Get amino acid compositions and Zc for genomes
  aa <- read.csv("MCK+23/genome_aa.csv")
  Zcvals <- canprot::Zc(aa)
  # List Zc for each genome in list
  Zclist <- lapply(genomes, function(genome) Zcvals[aa$organism %in% genome])
  # Use colors from Mateos et al., 2023
  dsr <- "#9c92ae" # "#bcb2ce"
  sox <- "#45b78d"
  mdd <- "#c24a96"
  # Colors for protein groups
  col <- c(dsr, sox, sox, dsr, mdd, mdd, mdd)

  # Plot Zc of genomes with S-cycling genes from Mateos et al. (2023)  20240916
  sulfur_Zc <- function() {
    # Setup plot
    # Need to make some adjustments after plotting with CHNOSZ::diagram()
    opar <- par(mar = c(4.1, 4.1, 2.1, 2.1), mgp = c(3.1, 1, 0), tcl = -0.5, yaxs = "r")
    n <- length(Zclist)
    boxplot(Zclist, col = col, names = character(n), xlab = "Age of earliest gene event (Ga)", ylab = "Zc of all proteins in genome", ylim = c(-0.25, -0.08))
    text(2.1, -0.24, CHNOSZ::hyphen.in.pdf("Sulfate-sulfite-sulfide"), col = dsr)
    text(2.8, -0.115, CHNOSZ::hyphen.in.pdf("Sulfate-\nthiosulfate"), col = sox)
    text(6.2, -0.22, "Organic sulfur", col = mdd)
    # Ages from Table 2 of Mateos et al.
    ages <- c("3.3-3.35  ", "  2.65-2.88", "2.6", "2.33-2.47", "2.28", "1.77", "0-2.38")
    # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
    text(x = 1:n, y = par()$usr[3] - 1.5 * strheight("A"), labels = ages, srt = 30, adj = 1, xpd = TRUE)
    genelab <- names(Zclist)
    genelab[genelab == "soxABXYZ"] <- "ABXYZ"
    axis(3, at = 1:n, labels = genelab, line = -0.8, lwd = 0, font = 3, gap.axis = -1)
    axis(3, at = 3, labels = "sox", line = 0.1, lwd = 0, font = 3)
    # Add number of genomes to labels
    n_genomes <- paste0("(", sapply(Zclist, length), ")")
    axis(3, at = 1:n, labels = n_genomes, line = -2, lwd = 0)
    #title(CHNOSZ::hyphen.in.pdf("Sulfur-cycling gene or gene cluster"), font.main = 1, line = 3)
    par(opar)
  }

  # Affinity ranking for genomes with different S-cycling genes
  sulfur_affinity <- function(panel) {
    if(is.null(panel)) par(mar = c(4.1, 4.1, 2.1, 4.1), mgp = c(2.5, 1, 0))
    basis("QEC+")
    # Keep genomes with single sulfur-cycling genes listed above
    myaa <- aa[aa$organism %in% unlist(genomes), ]
    # Load proteins for each genome
    ip <- add.protein(myaa, as.residue = TRUE)
    # Calculate affinity of composition reactions as a function of Eh and pH
    a <- affinity(pH = c(3, 10, res), O2 = c(-72.5, -58, res), iprotein = ip)
    # Group genomes according to presence of sulfur-cycling genes
    groups <- sapply(genomes, function(genome) match(genome, myaa$organism))
    # Calculate normalized sum of ranks for each group and make diagram
    amean <- agg.affinity(a, groups)
    # Lighten colors
    fill <- adjustcolor(col, alpha.f = 0.3)
    # Adjust labels
    names <- names(genomes)
    dx <- dy <- rep(0, length(names))
    dx[names == "soxABXYZ"] <- 1
    dy[names == "soxABXYZ"] <- -0.2
    # We need balance = 1 here to balance on residues 20250627
    diagram(amean, fill = fill, lty = 1, lwd = 2, font = 3, names = names, cex.names = 0.8, dx = dx, dy = dy, col = "gray20", balance = 1)
  }

  if("D" %in% panels) {
    sulfur_Zc()
    if(is.null(panel)) label.figure("D", font = 2, cex = 1.6)
  }
  if("E" %in% panels) {
    sulfur_affinity(panel)
    if(is.null(panel)) {
      # Add Eh7 axis
      stability_comparison(add = TRUE, pHlim = c(3, 10), Eh7_las = 0, datasets = numeric())
      title(main = CHNOSZ::hyphen.in.pdf("All proteins in genomes with specific S-cycling genes"), font.main = 1)
      label.figure("E", font = 2, cex = 1.6)
    } else {
      # Only add Eh7 axis
      stability_comparison(res = res, add = TRUE, Eh7_las = 0, datasets = numeric())
    }
  }

  if(pdf & is.null(panel)) dev.off()

}

# Comparison of Rubiscos and methanogen and Nitrososphaeria genomes
stability_comparison <- function(res = 400, add = FALSE, lwd = 2, lty = 1, pHlim = c(4, 10), O2lim = c(-72.5, -58), alpha.f = 1, Eh7_las = 1, datasets = 1:3) {

  # Setup basis species
  basis("QEC+")

  for(i in datasets) {

    if(i == 1) {
      # Methanogen proteomes 20220424
      # Read amino acid composition and compute Zc
      aa <- read.csv("DBCS23/methanogen_AA.csv")
      # Indices of Class I and Class II methanogens
      iI <- 20:36
      iII <- 1:19
      # Get the species in each group
      groups <- list("Class I" = iI, "Class II" = iII)
      col <- 7
      add <- add
    }

    if(i == 2) {
      # Thaumarchaeota (now Nitrososphaeria) 20220414
      # Amino acid compositions of predicted (Glimmer) and database (NCBI or IMG) proteomes
      predicted <- read.csv("DBCS23/Thaumarchaeota_predicted_AA.csv")
      database <- read.csv("DBCS23/Thaumarchaeota_database_AA.csv")
      # If both are available, use predicted instead of database
      aa <- rbind(predicted, database)
      aa <- aa[!duplicated(aa$organism), ]
      groupnames <- c("Basal", "Terrestrial", "Shallow", "Deep")
      # Get the species in each group
      groups <- sapply(groupnames, function(group) aa$protein == group, simplify = FALSE)
      # Compare Basal to Terrestrial 20240802
      groups <- groups[1:2]
      col <- 2
      add <- TRUE
    }

    if(i == 3) {
      # Rubisco 20240802
      # Read amino acid compositions
      fasta_file <- "KHAB17/rubisco.fasta"
      aa <- canprot::read_fasta(fasta_file)
      # Assign protein names
      aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)
      # Set up groups for affinity ranking:
      # 3 pre-GOE and 3 post-GOE proteins
      groups = list(pre = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), post = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
      col <- 4
      add <- TRUE
    }

    # Make affinity ranking plot 20220602
    # Load proteins and calculate affinity
    ip <- add.protein(aa, as.residue = TRUE)
    aout <- affinity(pH = c(pHlim, res), O2 = c(O2lim, res), iprotein = ip)
    # Calculate average ranking for each group and make diagram
    amean <- agg.affinity(aout, groups)
    diagram(amean, col = adjustcolor(col, alpha.f = alpha.f), lwd = lwd, lty = lty, add = add, names = "", balance = 1)

  }

  # Add Eh7 axis 20241218
  # Calculate range of pe at pH = 7
  # H2O = 0.5 O2(gas) + 2 H+ + 2 e- 
  # logK = -41.55
  # --> pe7 = 0.25 * logfO2 + 13.775
  # --> logfO2 = 4 * pe7 - 55.1
  Eh7_to_logfO2 <- function(Eh7) {
    pe7 <- convert(Eh7, "pe")
    logfO2 <- 4 * pe7 - 55.1
    logfO2
  }
  Eh7ticks <- seq(-0.25, 0.05, 0.05)
  logfO2ticks <- Eh7_to_logfO2(Eh7ticks)
  axis(4, at = logfO2ticks, labels = Eh7ticks, tcl = -0.3, mgp = c(2, 0.5, 0))
  mtext("Eh7 (V)", 4, line = 3, las = Eh7_las)

}

# Calculate average affinities for species in different groups
# 20220416 jmd first version (rank.affinity)
# 20250626 use aggregate function (mean or median) instead of average rank
agg.affinity <- function(aout, groups, fun = "mean") {

  # Put the affinities into matrix form
  amat <- sapply(aout$values, as.numeric)
  # Keep track of empty groups
  is_empty_group <- logical()

  # Get the average affinity for species in each group
  agg_values <- sapply(groups, function(group) {

    # Get number of species in this group
    if(inherits(group, "logical")) n <- sum(group)
    if(inherits(group, "integer")) n <- length(group)
    # Also handle indices classed as numeric 20250522
    if(inherits(group, "numeric")) n <- length(group)
    # Aggregate affinities
    group_values <- apply(amat[, group, drop = FALSE], 1, fun)

    # Remember empty group 20250527
    if(n == 0) {
      is_empty_group <<- c(is_empty_group, TRUE)
    } else {
      is_empty_group <<- c(is_empty_group, FALSE)
    }
    group_values

  })

  # Remove empty groups 20250527
  if(any(is_empty_group)) {
    agg_values <- agg_values[, !is_empty_group, drop = FALSE]
    empty_groups <- names(groups)[is_empty_group]
    message(paste("aggregate.affinity: removing empty groups:", paste(empty_groups, collapse = ", ")))
    groups <- groups[!is_empty_group, drop = FALSE]
  }

  # Restore dims
  dims <- dim(aout$values[[1]])
  if(getRversion() < "4.1.0") {
    # Using 'simplify = FALSE' in R < 4.1.0 caused error: 3 arguments passed to 'dim<-' which requires 2
    alist <- lapply(lapply(apply(agg_values, 2, list), "[[", 1), "dim<-", dims)
  } else {
    # apply() got 'simplify' argument in R 4.1.0 20230313
    alist <- apply(agg_values, 2, "dim<-", dims, simplify = FALSE)
  }
  aout$values <- alist

  # Rename species to group names (for use by diagram())
  aout$species <- aout$species[1:length(groups), ]
  aout$species$name <- names(groups)
  aout

}


##################
### SI Figures ###
##################

# Plot Zc of ancestral and extant nitrogenases from Garcia et al. (2020)  20250325
genoGOE_S1 <- function(ylim = c(-0.20, -0.12)) {

  ## Read FASTA file of ancient and extant sequences,
  ## downloaded from https://github.com/kacarlab/AncientNitrogenase.git
  #aa <- canprot::read_fasta("GMKK20/Extant-MLAnc_Align.fasta")
  # Get amino acid sequences precomputed from Extant-MLAnc_Align.fasta
  aa <- read.csv("GMKK20/nitrogenase_aa.csv")

  # List forms of nitrogenase and their ancestors
  form_to_anc <- list(
    "Clfx" = "D",
    "F-Mc" = "C",
    "Mb-Mc" = "B",
    "Anf" = "A",
    "Vnf" = "A",
    "Nif-II" = "E",
    "Nif-I" = "E"
  )

  # Start plot
  plot(extendrange(c(1, 7)), c(-0.20, -0.12), xlab = "Form of nitrogenase", xaxt = "n", ylab = "Zc", type = "n", ylim = c(-0.20, -0.12))
  axis(side = 1, at = seq_along(form_to_anc), labels = CHNOSZ::hyphen.in.pdf(names(form_to_anc)), gap.axis = 0)

  # Loop over nitrogenase forms
  set.seed(42)
  for(iform in seq_along(form_to_anc)) {

    # Calculate Zc of the ancestral proteins
    node <- form_to_anc[[iform]]
    ianc <- grepl(paste0("^Anc", node), aa$protein)
    Zc_anc <- canprot::Zc(aa[ianc, ])
    # Plot points with jitter
    xvals <- jitter(rep(iform, length(Zc_anc)), amount = 0.1)
    points(xvals, Zc_anc, pch = 19)

    # Calculate Zc of the extant proteins
    form <- names(form_to_anc[iform])
    iext <- which(aa$ref == form)
    Zc_ext <- canprot::Zc(aa[iext, ])
    xvals <- jitter(rep(iform, length(Zc_ext)), amount = 0.1)
    points(xvals, Zc_ext)

  }

  # Add legend and title
  legend("bottomright", c("Extant", "Ancestral"), pch = c(1, 19), bty = "n")

}
