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
  df <- read.csv("BacDive/cleaned_data_with_Zc.csv", check.names = FALSE)

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
    if(is.null(panel)) CHNOSZ::label.figure("A", font = 2, cex = 1.6)
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
    if(is.null(panel)) CHNOSZ::label.figure("B", font = 2, cex = 1.6)
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
      if(is.null(panel)) CHNOSZ::label.figure("C", font = 2, cex = 1.6, yfrac = 0.9)
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
    axis(3, at = c(1, 3, 6), labels = c("Entire genomes", "Binned by GC content", "Binned by metabolic cost"), tick = FALSE, font = 2)

    if(is.null(panel)) CHNOSZ::label.figure("D", font = 2, cex = 1.6, xfrac = 0.018)

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
    CHNOSZ::label.plot(LETTERS[j], font = 2, cex = 1.2, xfrac = 0.04)
  }
  # Outer axis labels
  mtext(ylab, side = 2, line = 3, adj = -0.68, font = 2)
  mtext("Gene\nage\n(Ma)", side = 4, line = 1.5, font = 2, las = 1, at = -0.0215, adj = 0.5)

  if(pdf) dev.off()

}

# Figure 5: Carbon oxidation state of reconstructed ancestral sequences and extant proteins 20250625
# New version of Figure 5 with stacked Zc plots above O2 profile  20260715 jmd
genoGOE_5 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure_5.pdf", width = 8, height = 8)
  # Setup figure region
  par(mfrow = c(6, 1))
  par(mar = c(0, 4.1, 0, 2.1))
  # Make plots
  #plot_phylostrata()
  plot_rubisco()
  plot_nitrogenase()
  plot_thioredoxin()
  plot_IPMDH()
  plot_oxygen()
  plot_temperature()
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
  legend("topright", lineages, pch = substr(lineages, 1, 1), col = coltext, bty = "n", cex = 1.1)
  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("Age (Ga)", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("Gene age groups", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("(Phylostrata: Liebeskind et al., 2016)", side = 3, line = -1.5, adj = 0.25, cex = par("cex") * 1.2)

}

# Plot Zc of rubisco from Kaçar et al. (2017)  20240407
plot_rubisco <- function() {
  # Read amino acid compositions
  fasta_file <- "KHAB17/rubisco.fasta"
  aa <- canprot::read_fasta(fasta_file)

  # Branch lengths from Fig. 4 of Kaçar et al. (2017) (root node at zero)
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
  labels <- round(predict(age2sub, newdata = data.frame(age = at)), 2)
  axis(1, at = at, labels = labels, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-0.20, -0.10, 0.05), labels = FALSE)
  axis(2, c(-0.20, -0.10), tick = FALSE, las = 1)

  # Add points and lines
  pch <- get_stages("rubisco", aa, return.pch = TRUE)
  points(ages, Zc_vals, pch = pch, bg = "black")
  lines(ages, Zc_vals, type = "b", pch = NA, col = 3)

  # Add text labels
  aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)
  text(ages[1:4], Zc_vals[1:4] + 0.007, aa$protein[1:4])
  text(ages[5], Zc_vals[5] + 0.007, aa$protein[5], adj = 1)
  text(ages[6], Zc_vals[6] - 0.007, aa$protein[6], adj = 0)

  # Add legend
  legend("topright", legend = "Linear age model\nApproximate ages only", bty = "n", text.font = 3, cex = 1.1)
  legend("bottomright", legend = paste("Stage", 1:3), pch = c(21, 22, 24), pt.bg = "black", title = "Stages for\nthermodynamic analysis", bty = "n")
  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("aa subs/site", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("A. Rubisco large subunit", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  legend("topleft", "RAS: Kaçar et al. (2017)", bty = "n", title = "", inset = c(-0.023, 0), cex = 1.1)
}

# Plot Zc for ancestral and modern nitrogenases from Cuevas Zuviría et al. (2025)  20260715
# Add data from Rucker et al. (2026)  20260801
plot_nitrogenase <- function() {

  # List of all proteins with branch lengths
  # (aa_sub - measured from root node in CDA+25 Fig. 3)
  # * Anc1209 and Anc1224 are excluded because they are similar to Anc4 and Anc3 in RBH+26
  # N | CDA+25   | RBH+26    | aa_sub
  # 1 | Anc1206  |           | 0.22
  # 2 | Anc1207  |           | 0.71
  #   | Anc1209* |           | 0.87
  # 3 |          | Anc4_1223 | 0.93
  #   | Anc1224* |           | 1.04
  # 4 |          | Anc3_1231 | 1.12
  # 5 |          | Anc2_1304 | 1.48
  # 6 |          | Anc1_1312 | 1.64
  # 7 | A. vinelandii        | 1.76

  # Node numbers for both datasets (lower = older)
  x_CDA <- c(1, 2, 7)
  x_RBH <- c(3, 4, 5, 6, 7)

  # CDA+25 node names
  CDA_nodes <- c("1206_map", "1207_map", "Nif_Azotobacter_vinelandii")
  # RBH+26 node names
  RBH_nodes <- c("Anc4_1223", "Anc3_1231", "Anc2_1304", "Anc1_1312", "WT_Nif_Azotobacter_vinelandii")
  # Combined node names from both datasets
  all_nodes <- character()
  all_nodes[x_CDA] <- CDA_nodes
  all_nodes[x_RBH] <- RBH_nodes

  # Branch lengths for all proteins
  aa_sub <- c(0.22, 0.71, 0.93, 1.12, 1.48, 1.64, 1.76)
  # Ages for proteins (Anc4 - Anc1) in RBH+26 Fig. 1
  age_RBH <- c(2.00, 1.72, 1.35, 1.17)

  # Linear fit between branch lengths and ages from RBH+26
  aa_sub_RBH <- aa_sub[head(x_RBH, -1)]
  sub2age <- lm(age_RBH ~ aa_sub_RBH)
  print(paste("R-squared for age_RBH ~ aa_sub_RBH:", round(summary(sub2age)$r.squared, 3)))
  # Predict ages for all proteins from branch length
  ages <- predict(sub2age, newdata=data.frame(aa_sub_RBH = aa_sub))
  # Make inverse model to get tick labels
  age2sub <- lm(aa_sub_RBH ~ age_RBH)

  # Start plot
  ylim <- c(-0.20, -0.12)
  xlim <- c(4.5, 0)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  # Put labels at whole Ga values
  at <- c(3, 2, 1)
  labels <- round(predict(age2sub, newdata = data.frame(age_RBH = at)), 1)
  axis(1, at = at, labels = labels, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-0.20, -0.12, 0.04), labels = FALSE)
  axis(2, c(-0.20, -0.12), tick = FALSE, las = 1)

  # Loop over Nif subunits
  subunits <- c("D", "K", "H")
  for(i in 1:length(subunits)) {
    # Initialize Zc values
    Zc_vals <- numeric()
    # Get CDA+25 data
    fasta_file <- paste("AGNifAlign105.ext-anc.alt", subunits[i], "fasta", sep = ".")
    fasta_path <- file.path("CDA+25", fasta_file)
    aa_CDA <- canprot::read_fasta(fasta_path)
    # Get RBH+26 data
    fasta_file <- paste0("Nif", subunits[i], "_selected_seqs.fasta")
    fasta_path <- file.path("RBH+26", fasta_file)
    aa_RBH <- canprot::read_fasta(fasta_path)
    # Combine the data frames
    aa <- rbind(aa_CDA, aa_RBH)
    # Find the rows with matching node names
    iaa <- match(all_nodes, aa$protein)
    node_aa <- aa[iaa, ]
    # Calculate Zc
    Zc_vals <- canprot::Zc(node_aa)
    # Add points and lines
    pch <- get_stages("nitrogenase", node_aa, return.pch = TRUE)
    points(ages, Zc_vals, pch = pch, bg = "black")
    lines(ages, Zc_vals, type = "b", pch = NA, col = 4)
    # Add subunit labels
    text(ages[1], Zc_vals[1], subunits[i], adj = 1.8)
    text(ages[7], Zc_vals[7], subunits[i], adj = -0.8)
  }

  # Add legend
  legend("topright", legend = "Linear age model\nApproximate ages only", bty = "n", text.font = 3, cex = 1.1)
  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("aa subs/site", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("B. Nitrogenase subunits", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  legend("topleft", "RAS: Cuevas Zuviría et al. (2025)\nand Rucker et al. (2026)", bty = "n", title = "", inset = c(-0.023, 0), cex = 1.1)
}

# Plot Zc for ancestral thioredoxins from Perez-Jimenez et al. (2011)  20250625
plot_thioredoxin <- function() {
  # Read data file with ages and PDB IDs from Del Galdo et al. (2019)
  dat <- read.csv("PIZ+11/DAAD19.csv")
  # Read amino acid compositions
  aa <- canprot::read_fasta("PIZ+11/thioredoxin.fasta")
  aa$protein <- dat$name
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
  for(lineage in c("Bacteria", "Arc-Euk")) {
    ilineage <- dat$Lineage == lineage
    # Add points and lines
    pch <- get_stages("thioredoxin", aa[ilineage, ], return.pch = TRUE)
    points(dat$Age[ilineage], Zc[ilineage], pch = pch, bg = "black")
    lines(dat$Age[ilineage], Zc[ilineage], type = "b", pch = NA, col = 7)
    add_age_error_bars(ilineage, Zc)
  }
  text(3.5, -0.22, "Bacteria")
  text(2.5, -0.26, "Archaea+Eukaryota")

  # Add axis labels and title
  mtext("Zc", side = 2, las = 1, line = 1.5, font = 2, cex = par("cex") * 1.2)
  mtext("Age (Ga)", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("C. Thioredoxin", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  legend("topleft", CHNOSZ::hyphen.in.pdf("RAS: Perez-Jimenez et al. (2011)"), bty = "n", title = "", inset = c(-0.022, 0), cex = 1.1)
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

  # Add points and lines
  pch <- get_stages("IPMDH", aa, return.pch = TRUE)
  points(ages / 1000, Zc, pch = pch, bg = "black")
  lines(ages / 1000, Zc, type = "b", pch = NA, col = 5)

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
  mtext(CHNOSZ::hyphen.in.pdf("D. 3-isopropylmalate dehydrogenase (IPMDH)"), side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  legend("topleft", "RAS: Cui et al. (2025)", bty = "n", title = "", inset = c(-0.025, 0), cex = 1.1)
}

# Plot O2 curve from Lyons et al. (2024)  20260716
plot_oxygen <- function() {
  # Ribbon color
  col <- "#3E7CB1"

  # Start plot
  xlim <- c(4.5, 0)
  ylim <- c(-6, 0.2)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  # Axis ticks and labels
  axis(1, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(-4, 0, 1), labels = FALSE)
  axis(2, c(-4, -2, 0), tick = FALSE, las = 1)
  # Ticks below breaks in scale
  axis(2, c(-5, -5.75), labels = c(-7, -13), las = 1)
  # Scale break marks
  text(par("usr")[1], -4.6, "~", xpd = NA, cex = 1.5, srt = 15)
  text(par("usr")[1], -5.5, "~", xpd = NA, cex = 1.5, srt = 15)

  # GOE, NOE, and Devonian
  rect(2.4, -5, 2, par("usr")[4], border = NA, col = "#d0edfd")
  text(2.2, -5.55, "GOE", cex = 1.1)
  rect(0.8, -5, 0.54, par("usr")[4], border = NA, col = "#d0edfd")
  text(0.67, -5.55, "NOE", cex = 1.1)
  rect(0.419, -5, 0.359, par("usr")[4], border = NA, col = "#d0edfd")
  text(0.389, -5.55, "  Devonian", cex = 1.1)
  box()

  ## Archean lines
  #lines(c(4, 3), c(-5.6, -5.6), col = col, lwd = 15)
  #lines(c(3, 2.45), c(-5, -5), col = col, lwd = 15)

  # Archean lines, with a cross-fade between the two bars centered at 3.0 Ga
  fade_range <- c(3.5, 2.5)
  fade_offset <- 0.3
  fade_x <- seq(fade_range[1], fade_range[2], length.out = 100)
  lines(c(4, fade_range[1]), c(-5.6, -5.6), col = col, lwd = 15)
  lines(c(fade_range[2], 2.45), c(-5, -5), col = col, lwd = 15)
  for (i in seq_len(length(fade_x) - 1)) {
    frac <- (i - 1) / (length(fade_x) - 1)
    if(i > 25) lines(fade_x[c(i, i + 1)] + fade_offset, c(-5.6, -5.6), col = adjustcolor(col, alpha.f = 1 - frac), lwd = 15)
    if(i < 75) lines(fade_x[c(i, i + 1)] - fade_offset, c(-5, -5), col = adjustcolor(col, alpha.f = frac), lwd = 15)
  }

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
  mtext("E. Atmospheric oxygen", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  legend("topleft", "Adapted from Lyons et al. (2024)", bty = "n", title = "", inset = c(-0.0215, 0), cex = 1.1)
}

# Plot temperature curves from IR24 and JSW07  20260727 jmd
plot_temperature <- function() {
  # Start plot
  xlim <- c(4.5, 0)
  ylim = c(0, 40)
  plot(xlim, ylim, xlim = xlim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", type = "n")
  # Axis ticks and labels
  axis(1, tcl = 0.5, mgp = c(-3, -1.5, 0))
  axis(2, seq(0, 40, 10), labels = FALSE)
  axis(2, c(0, 40), tick = FALSE, las = 1)

  # Add line from IR24
  dat <- read.csv("IR24/IR24_Fig3B.csv")
  lines(dat$Age_Ga, dat$T_C, col = 4, lwd = 2)

  # Add line from JSW07
  dat <- read.csv("JSW07/JSW07_Fig15.csv")
  lines(dat$Age_Ma/1000, dat$T_C, col = 2, lwd = 1.5, lty = 2)

  # Identify min/max values
  T_range <- range(dat$T_C)
  Age_vals <- dat$Age_Ma[dat$T_C %in% T_range]/1000
  print(paste("Age range for min/max T in IR24:", paste(round(Age_vals, 2), collapse = ", ")))
  T_text <- sapply(T_range, function(T_val) {
    bquote(.(round(T_val)) ~ degree*C)
  })
  text(Age_vals, T_range, as.expression(T_text), adj = -0.2)

  # Add loess curve for JSW07
  dat.lo <- loess(T_C ~ Age_Ma, dat)
  lines(dat$Age_Ma/1000, predict(dat.lo), lty = 3)

  # Add axis labels and title
  mtext(quote(bolditalic(T)~bold("("*degree*C*")")), side = 2, las = 1, line = 1.2, font = 2, cex = par("cex") * 1.2)
  mtext("Age (Ga)", side = 1, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)
  mtext("F. Ocean surface temperature", side = 3, line = -1.5, font = 2, adj = 0.01, cex = par("cex") * 1.2)

  # Add legend
  legend("topleft", c("Jaffrés et al. (2007)", "Loess fit", "Isson and Rauzi (2024)"),
    lty = c(2, 3, 1), col = c(2, 1, 4), lwd = c(1.5, 1.5, 1.5), title = "", bty = "n", cex = 1.1)
}


#####################################
### End of functions for Figure 5 ###
#####################################

# Figure 6: Relative stabilities of proteins as a function of environmental variables
genoGOE_6 <- function(pdf = FALSE, panel = NULL) {

  if(is.null(panel)) {
    if(pdf) pdf("Figure_6.pdf", width = 11, height = 9)
    mat <- matrix(c(1,1,1,  2,2,2, 3,
                    4,4, 5,5, 6,6, 7,
                    8,8,8,  9,9,9, 10),
                  nrow = 3, byrow = TRUE)
    layout(mat, heights = c(5, 4, 5), widths = c(1,1,1, 1,1,1, 0.5))
    par(cex = 1)
  }
  panels <- if(is.null(panel)) LETTERS[1:5] else panel
  # Margin setting for CHNOSZ::diagram() (via plot_stability())
  mar <- c(3, 3.5, 2, 1)

  # Panel A: Compare Zc of ancestral Rubisco sequences from different sources  20260717
  if("A" %in% panels) {
    # Names of ancestral proteins
    anc_names <- list(
      KHAB17 = c("Anc_I/II/III", "Anc_I/III", "Anc_I/III'", "Anc_I", "Anc_IA/B", "Anc_IB"),
      "SGZ+22" = c("AncL", "AncLS"),
      ACK25 = c("Anc_I_III", "Anc_I_I-prime", "Anc_I-prime", "Anc_I", "Anc_IAB", "Anc_IB", "Anc_ICD", "Anc_IA")
    )
    # Relative ages of ancestral proteins (lower = older)
    anc_ages <- list(
      KHAB17 = c(1, 2, 3, 6, 7, 8),
      "SGZ+22" = c(4, 6),
      ACK25 = c(2, 4, 5, 6, 7, 8, 9, 10)
    )
    # Start plot
    xlab <- "Ancestral nodes (older to younger)"
    xlim <- c(1, 10)
    ylim <- c(-0.20, -0.10)
    opar <- par(mar = c(4, 4, 2, 1))
    plot(xlim, ylim, type = "n", xlab = xlab, ylab = "Zc", ylim = ylim, xaxt = "n", las = 1)
    # Add GOE range (Model 1) 20260723
    rect(4, par("usr")[3], 6, par("usr")[4], border = NA, col = "#d0edfd")
    text(5, -0.19, "GOE\n(Model 1)")
    box()
    # Add tick marks
    all_nodes <- c("I/II/III", "I/III", "I/III'", "I/I'\nAncL", "I'", "I\nAncLS", "IAB", "IB", "ICD", "IA")
    axis(1, 1:10, all_nodes, padj = 1, mgp = c(3, 0, 0), gap.axis = 0)

    # Plot Zc of rubisco from Kaçar et al. (2017)  20240407
    # Read amino acid compositions
    fasta_file <- "KHAB17/rubisco.fasta"
    aa <- canprot::read_fasta(fasta_file)
    # Make sure proteins are ordered as expected
    stopifnot(all.equal(aa$protein, anc_names$KHAB17))
    # Get ages
    ages <- anc_ages$KHAB17
    # Get Zc values
    Zc_vals <- canprot::Zc(aa)
    # Add points
    lines(ages, Zc_vals, lwd = 2, cex = 1.5, pch = 19, col = 2, type = "b")
    # Plot Zc for ancestral sequences from Amritkar et al., 2025  20260713
    seq_dir <- "ACK25"
    # Get Zc of LSU ancestors
    Zc_LSU <- numeric()
    for(LSU_Anc in anc_names$ACK25) {
      seq_file <- paste0(file.path(seq_dir, LSU_Anc), "_LSU.fasta")
      # Use read_fasta() to get amino acid composition
      aa <- canprot::read_fasta(seq_file)
      # Calculate and save Zc
      Zc_LSU <- c(Zc_LSU, canprot::Zc(aa))
    }
    # Add points
    ages <- anc_ages$ACK25
    lines(ages, Zc_LSU, lwd = 2, cex = 1.5, type = "b")
    # Plot Zc for ancestral sequences from Schulz et al., 2022  20260713
    aa <- canprot::read_fasta("SGZ+22/sequences.fasta")
    stopifnot(all.equal(aa$protein, anc_names$"SGZ+22"))
    Zc_vals <- canprot::Zc(aa)
    lines(anc_ages$"SGZ+22", Zc_vals, lwd = 2, cex = 1.5, pch = 15, col = 4, type = "b")
    # Add lines and text for stages
    abline(v = c(3.5, 6.5), lty = 2, lwd = 2, col = 8)
    text(2, -0.09, "Stage 1", font = 3, col = 8, xpd = NA)
    text(5, -0.09, "Stage 2", font = 3, col = 8, xpd = NA)
    text(8.5, -0.09, "Stage 3", font = 3, col = 8, xpd = NA)
    # Add GOE Model 2 text
    text(6.6, -0.19, "Start GOE\n(Model 2)", adj = 0)
    # Add legend
    legend("topleft", c("Amritkar", "Schulz", "Kaçar"), pch = c(1, 15, 19), col = c(1, 4, 2), pt.cex = 1.5, bty = "n")
    CHNOSZ::label.figure("A", cex = 1.5, font = 2, yfrac = 0.94)
    par(opar)
  }

  if("B" %in% panels) {
    # Rubisco relative stability
    plot_stability("rubisco_1_2", plot_names = FALSE, mar = mar, col = 3)
    plot_stability("rubisco_2_3", plot_names = FALSE, add = TRUE, col = 3)
    text(7, -64.5, "Stage 1")
    text(7, -62.7, "Stage 2")
    text(6, -68.5, "Stage 2")
    text(6, -66.3, "Stage 3")
    title("Rubisco", font.main = 1)
    CHNOSZ::label.figure("B", cex = 1.5, font = 2, yfrac = 0.94)
  }

  # Arrow for Rubisco redox evolution
  if(is.null(panel)) {
    opar <- par(mar = c(3, 0, 2, 0))
    plot.new()
    plot.window(c(0, 1), c(-72.5, -58))
    arrows(0.0, -63.2, 0, -70, length = 0.15, xpd = NA)
    text(0.1, -66.6, "Time", adj = 0, font = 3)
    text(-0.1, -60.5, "Retrograde\nredox\nevolution", adj = 0, font = 3, xpd = NA)
    par(opar)
  }

  if("C" %in% panels) {
    arrowfun <- function(y1, y2) {
      arrows(9.5, y1, 9.5, y2, length = 0.15)
      text(9.3, y2, "Time", adj = 1, font = 3)
    }

    plot_stability("nitrogenase", mar = mar, col = 4)
    title("Nitrogenase (D and K subunits)", font.main = 1)
    CHNOSZ::label.figure("C", cex = 1.5, font = 2, yfrac = 0.94)

    plot_stability("thioredoxin", mar = mar, col = 7)
    title("Thioredoxin (all lineages)", font.main = 1)

    plot_stability("IPMDH", mar = mar, col = 5)
    title("IPMDH", font.main = 1)
  }

  # Arrow for other proteins' redox evolution
  if(is.null(panel)) {
    opar <- par(mar = c(3, 0, 2, 0))
    plot.new()
    plot.window(c(0, 1), c(-72.5, -58))
    arrows(0.0, -70, 0, -64, length = 0.15, xpd = NA)
    text(0.1, -67, "Time", adj = 0, font = 3)
    text(-0.1, -60.5, "Prograde\nredox\nevolution", adj = 0, font = 3, xpd = NA)
    par(opar)
  }

  # Evolutionary oxidation and relative stabilities for genomes with S-cycling genes 20241211
  # genoGOE/sulfur_genomes.R
  # 20241211 add Eh-pH affinity aggregation
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

  # Mean affinities for genomes with different S-cycling genes
  sulfur_affinity <- function(panel) {
    if(is.null(panel)) par(mar = c(4.1, 4.1, 2.1, 4.1), mgp = c(2.5, 1, 0))
    basis("QEC+")
    # Keep genomes with single sulfur-cycling genes listed above
    myaa <- aa[aa$organism %in% unlist(genomes), ]
    # Load proteins for each genome
    ip <- add.protein(myaa, as.residue = TRUE)
    # Calculate affinity of composition reactions as a function of Eh and pH
    res <- 300
    a <- affinity(pH = c(3, 10, res), O2 = c(-72, -62, res), iprotein = ip)
    # Group genomes according to presence of sulfur-cycling genes
    groups <- sapply(genomes, function(genome) match(genome, myaa$organism))
    # Calculate mean affinities for each group and make diagram
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
    if(is.null(panel)) CHNOSZ::label.figure("D", font = 2, cex = 1.6, yfrac = 0.94)
  }
  if("E" %in% panels) {
    sulfur_affinity(panel)
    # Add Eh7 axis
    add_Eh7_axis(las = 0)
    title(main = CHNOSZ::hyphen.in.pdf("        All proteins in genomes with specific S-cycling genes"), font.main = 1)
    CHNOSZ::label.figure("E", font = 2, cex = 1.6, yfrac = 0.94)
  }

  if(pdf & is.null(panel)) dev.off()

}

# Relative stability diagrams for reconstructed ancestral proteins 20260721
plot_stability <- function(dataset = "rubisco", res = 200, pHlim = c(4, 10), O2lim = c(-72, -62), add = FALSE, plot_names = TRUE, mar = NULL, col = 1) {

  # Setup basis species
  basis("QEC+")

  if(dataset == "methanogen_old") {
    # Methanogen proteomes 20220424
    # Read amino acid composition and compute Zc
    aa <- read.csv("DBCS23/methanogen_AA.csv")
    # Indices of Class I and Class II methanogens
    iI <- 20:36
    iII <- 1:19
    # Get the species in each group
    stages <- list("Class I" = iI, "Class II" = iII)
  }

  if(dataset == "thaumarchaeota_old") {
    # Thaumarchaeota (now Nitrososphaeria) 20220414
    # Amino acid compositions of predicted (Glimmer) and database (NCBI or IMG) proteomes
    predicted <- read.csv("DBCS23/Thaumarchaeota_predicted_AA.csv")
    database <- read.csv("DBCS23/Thaumarchaeota_database_AA.csv")
    # If both are available, use predicted instead of database
    aa <- rbind(predicted, database)
    aa <- aa[!duplicated(aa$organism), ]
    stagenames <- c("Basal", "Terrestrial", "Shallow", "Deep")
    # Get the species in each group
    stages <- sapply(stagenames, function(stage) aa$protein == stage, simplify = FALSE)
    # Compare Basal to Terrestrial 20240802
    stages <- stages[1:2]
  }

  if(dataset == "rubisco_old") {
    # Rubisco 20240802
    # Read amino acid compositions
    fasta_file <- "KHAB17/rubisco.fasta"
    aa <- canprot::read_fasta(fasta_file)
    # Assign protein names
    aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)
    # Set up groups for affinity aggregation:
    # 3 Form I/(II)/III and 3 Form I proteins (possibly pre-GOE and post-GOE)
    stages = list(
      pre = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
      post = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
    )
  }

  if(dataset %in% c("rubisco", "rubisco_1_2", "rubisco_2_3", "rubisco_KHAB17_1_2", "rubisco_KHAB17_2_3", "rubisco_ACK25")) {

    # Stages of Rubisco evolution using multiple datasets 20260721
    # Protein names in FASTA files
    anc_names <- list(
      KHAB17 = c("Anc_I/II/III", "Anc_I/III", "Anc_I/III'", "Anc_I", "Anc_IA/B", "Anc_IB"),
      "SGZ+22" = c("AncL", "AncLS"),
      ACK25 = c("Anc_I_III", "Anc_I_I-prime", "Anc_I-prime", "Anc_I", "Anc_IAB", "Anc_IB", "Anc_ICD", "Anc_IA")
    )
    # Evolutionary stages
    anc_ages <- list(
      KHAB17 = c(1, 1, 1, 2, 3, 3),
      "SGZ+22" = c(2, 2),
      ACK25 = c(1, 2, 2, 2, 3, 3, 3, 3)
    )

    # Dataset 1: Kacar et al., 2017
    aa1 <- canprot::read_fasta("KHAB17/rubisco.fasta")
    stopifnot(all.equal(aa1$protein, anc_names$KHAB17))
    # Dataset 2: Schulz et al., 2022
    aa2 <- canprot::read_fasta("SGZ+22/sequences.fasta")
    stopifnot(all.equal(aa2$protein, anc_names$"SGZ+22"))
    # Dataset 3: Amritkar et al., 2025
    aa3_list <- lapply(anc_names$ACK25, function(LSU_Anc){
      seq_file <- paste0(file.path("ACK25", LSU_Anc), "_LSU.fasta")
      canprot::read_fasta(seq_file)
    })
    aa3 <- do.call(rbind, aa3_list)
    aa3$protein <- gsub("_LSU", "", aa3$organism)
    aa3$organism <- "LSU"
    stopifnot(all.equal(aa3$protein, anc_names$ACK25))

    # Select single source or put together all sources
    if(grepl("KHAB17", dataset)) aa <- aa1
    else if(grepl("ACK25", dataset)) aa <- aa3
    else aa <- rbind(aa1, aa2, aa3)
    # Check that all protein-organism pairs are unique 
    # (so that none are dropped by add.protein())
    stopifnot(!any(duplicated(paste(aa$protein, aa$organism, sep = "_"))))

    # Set up groups for affinity aggregation
    stages <- get_stages("rubisco", aa)
    if(grepl("1_2", dataset)) stages <- stages[1:2]
    if(grepl("2_3", dataset)) stages <- stages[2:3]

  }

  if(dataset == "IPMDH") {
    # IPMDH 20260721
    aa <- canprot::read_fasta("CDY+25/IPMDH.fasta")
    stages <- get_stages("IPMDH", aa)
  }

  if(dataset %in% c("nitrogenase", "nitrogenase_1_2", "nitrogenase_2_3")) {
    # Get data from Cuevas Zuviría et al. (2025)
    aa_D_CDA <- canprot::read_fasta(file.path("CDA+25", "AGNifAlign105.ext-anc.alt.D.fasta"))
    aa_K_CDA <- canprot::read_fasta(file.path("CDA+25", "AGNifAlign105.ext-anc.alt.K.fasta"))
    # Get data from Rucker et al. (2026)
    aa_D_RBH <- canprot::read_fasta(file.path("RBH+26", "NifD_selected_seqs.fasta"))
    aa_K_RBH <- canprot::read_fasta(file.path("RBH+26", "NifK_selected_seqs.fasta"))
    aa <- rbind(aa_D_CDA, aa_K_CDA, aa_D_RBH, aa_K_RBH)

    # Set up groups for affinity aggregation
    stages <- get_stages("nitrogenase", aa)
    if(grepl("1_2", dataset)) stages <- stages[1:2]
    if(grepl("2_3", dataset)) stages <- stages[2:3]
  }

  if(dataset %in% c("thioredoxin", "thioredoxin_A", "thioredoxin_B", "thioredoxin_B_1_2", "thioredoxin_B_2_3")) {
    # Read data file with ages and PDB IDs from Del Galdo et al. (2019)
    dat <- read.csv("PIZ+11/DAAD19.csv")
    # Read amino acid compositions
    aa <- canprot::read_fasta("PIZ+11/thioredoxin.fasta")
    stopifnot(all.equal(sapply(strsplit(aa$protein, "_"), "[", 1), dat$PDB))
    aa$protein <- dat$name
    # Subset Arc-Euk and Bacteria
    if(dataset == "thioredoxin_A") aa <- aa[dat$Lineage == "Arc-Euk", ]
    if(dataset == "thioredoxin_B") aa <- aa[dat$Lineage == "Bacteria", ]
    stages <- get_stages("thioredoxin", aa)
    if(dataset == "thioredoxin_B_1_2") stages <- stages[1:2]
    if(dataset == "thioredoxin_B_2_3") stages <- stages[2:3]
  }

  # Make affinity plot 20220602
  # Load proteins and calculate affinity
  ip <- add.protein(aa, as.residue = TRUE)
  aout <- affinity(pH = c(pHlim, res), O2 = c(O2lim, res), iprotein = ip)
  # Calculate average affinity for each stage and make diagram
  amean <- agg.affinity(aout, groups = stages)
  lwd <- 2
  lty <- 1
  if(plot_names) names <- amean$species$name else names <- NA
  diagram(amean, col = col, lwd = lwd, lty = lty, add = add, balance = 1, names = names, format.names = FALSE, mar = mar)

}

# Function to get evolutionary stages for Fig. 5 and 6  20260723
get_stages <- function(dataset, aa, return.pch = FALSE) {
  if(dataset == "rubisco") {
    stages = list(
      "Stage 1" = aa$protein %in% c("Anc_I/II/III", "Anc_I/III", "Anc_I/III'", "Anc_I_III"),
      "Stage 2" = aa$protein %in% c("Anc_I", "AncL", "AncLS", "Anc_I_I-prime", "Anc_I-prime", "Anc_I"),
      "Stage 3" = aa$protein %in% c("Anc_IA/B", "Anc_IB", "Anc_IAB", "Anc_IB", "Anc_ICD", "Anc_IA")
    )
  } else if(dataset == "nitrogenase") {
    stages = list(
      "Stage 1" = aa$protein %in% c("1206_map", "1207_map"),
      "Stage 2" = aa$protein %in% c("Anc4_1223", "Anc3_1231"),
      "Stage 3" = aa$protein %in% c("Anc2_1304", "Anc1_1312", "Nif_Azotobacter_vinelandii", "WT_Nif_Azotobacter_vinelandii")
    )
  } else if(dataset == "IPMDH") {
    stages = list(
      "Stage 1" = aa$protein %in% c("Anc01", "Anc02", "Anc03", "Anc04"),
      "Stage 2" = aa$protein %in% c("Anc05", "Anc06", "Anc07", "Anc08"),
      "Stage 3" = aa$protein %in% c("Anc09", "Anc10", "Anc11", "EcIPMDH")
    )
  } else if(dataset == "thioredoxin") {
    stages = list(
      "Stage 1" = aa$protein %in% c("LBCA", "AECA", "LACA", "LPBCA"),
      "Stage 2" = aa$protein %in% c("LGPCA", "LECA", "LAFCA"),
      "Stage 3" = aa$protein %in% c("Ecoli", "Human")
    )
  }
  if(!return.pch) return(stages) else {
    # Return pch values for Fig. 5
    pch.ind <- rep(NA, nrow(aa))
    for(i in 1:length(stages)) pch.ind[stages[[i]]] <- i
    pch.vals <- c(21, 22, 24, 23)
    return(pch.vals[pch.ind])
  }
}

add_Eh7_axis <- function(las = 1) {
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
  mtext("Eh7 (V)", 4, line = 3, las = las, cex = par("cex"))
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

# Figure S1: Zc of ancestral and extant nitrogenases from Garcia et al. (2020)  20250325
genoGOE_S1 <- function(pdf = FALSE) {

  ## From scratch: Read FASTA file of ancient and extant sequences,
  ## available at https://github.com/kacarlab/AncientNitrogenase.git
  #aa <- canprot::read_fasta("GMKK20/Extant-MLAnc_Align.fasta")

  # Get amino acid sequences precomputed from Extant-MLAnc_Align.fasta
  aa <- read.csv("GMKK20/nitrogenase_aa.csv")
  # Use only proteins shown in Fig. 6 of Garcia et al. (2020)
  from_GMKK20_Fig6 <- sapply(aa$abbrv == "GMKK20_Fig6", isTRUE)
  aa <- aa[from_GMKK20_Fig6, ]

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
  if(pdf) pdf("Figure_S1.pdf", width = 6, height = 6)
  layout(matrix(c(1, 2, 1, 3), nrow = 2))
  par(mar = c(4, 4, 1, 1))
  ylim <- c(-0.20, -0.12)
  plot(extendrange(c(1, 7)), ylim, xlab = "Form of nitrogenase", xaxt = "n", ylab = "Zc", type = "n", ylim = ylim)
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
    aa_ext <- aa[iext, ]
    Zc_ext <- canprot::Zc(aa_ext)
    xvals <- jitter(rep(iform, length(Zc_ext)), amount = 0.1)
    points(xvals, Zc_ext)
  }
  # Add legend
  legend("bottomright", c("Extant", "Ancestral"), pch = c(1, 19), bty = "n")
  label.figure("A", cex = 1.5, font = 2, xfrac = 0.02)

  # Panel B: Cyanobacteria vs others for Nif-I 20260726
  par(mar = c(4, 4, 3, 1))
  # Read CSV again to get all proteins
  aa <- read.csv("GMKK20/nitrogenase_aa.csv")
  # Limit to Nif-I, treating NA as FALSE
  aa <- aa[sapply(aa$ref == "Nif-I", isTRUE), ]
  # Identify Cyanobacteria
  icyano <- sapply(aa$organism == "Cyanobacteriota", isTRUE)
  # Get Zc values
  Zc_list <- list(
    Cyano = canprot::Zc(aa[icyano, ]),
    Other = canprot::Zc(aa[!icyano, ])
  )
  names(Zc_list)[1] <- ""
  names(Zc_list)[2] <- ""
  boxplot(Zc_list, ylab = "Zc")
  # Plot x-axis labels with custom spacing
  names(Zc_list)[1] <- paste0("Cyanobacteriota\n(", length(Zc_list[[1]]), ")")
  names(Zc_list)[2] <- paste0("Other phyla\n(", length(Zc_list[[2]]), ")")
  axis(1, at = 1:2, labels = names(Zc_list), mgp = c(3, 2, 0))
  pval <- t.test(Zc_list[[1]], Zc_list[[2]], alternative = "greater")$p.value
  legend("topleft", legend = bquote(italic(p) == .(signif(pval, 2))), bty = "n")
  title(CHNOSZ::hyphen.in.pdf("Extant Nif-I"), font.main = 1)
  label.figure("B", cex = 1.5, font = 2, yfrac = 0.92)

  # Panel C: Nif-I vs Nif-II 20260726
  aa <- read.csv("GMKK20/nitrogenase_aa.csv")
  aa <- aa[sapply(aa$ref %in% c("Nif-I", "Nif-II"), isTRUE), ]
  inifI <- aa$ref == "Nif-I"
  Zc_list <- list(
    NifI = canprot::Zc(aa[inifI, ]),
    NifII = canprot::Zc(aa[!inifI, ])
  )
  names(Zc_list)[1] <- CHNOSZ::hyphen.in.pdf(paste0("Nif-I (", length(Zc_list[[1]]), ")"))
  names(Zc_list)[2] <- CHNOSZ::hyphen.in.pdf(paste0("Nif-II (", length(Zc_list[[2]]), ")"))
  boxplot(Zc_list, ylab = "Zc")
  pval <- t.test(Zc_list[[1]], Zc_list[[2]], alternative = "greater")$p.value
  legend("topright", legend = bquote(italic(p) == .(signif(pval, 2))), bty = "n")
  title(CHNOSZ::hyphen.in.pdf("Extant Nif-I vs Nif-II"), font.main = 1)
  label.figure("C", cex = 1.5, font = 2, yfrac = 0.92)

  if(pdf) dev.off()

}

# Figure S2: Relative stabilities for Rubisco (pairwise and groupwise affinities)
genoGOE_S2 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_S2.pdf", width = 7, height = 5)
  mat <- matrix(c(1,1, 2,2, 0, 3,3, 0), nrow = 2, byrow = TRUE)
  layout(mat)

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
  res <- 200
  
  # Panel A: Pairwise stability boundaries for Rubisco
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
        text(x, y - 0.04, aa$protein[pre], cex = 0.6)
        text(x, y + 0.04, aa$protein[post], cex = 0.6)
      }
    }
  }

  text(5.2, 0.62, "Higher affinity\nfor Form I protein\nin each pair")
  text(4.8, -0.15, "Higher affinity for\nForm I/(II)/III protein in each pair", srt = -21)
  title("Pairwise Rubiscos", font.main = 1)
  label.figure("A", cex = 1.5, font = 2, yfrac = 0.92)

  # Panel B: Groupwise stability boundary (Eh-pH diagram)
  # Calculate affinity of composition reactions for all proteins
  aout <- affinity(pH = c(0, 14, res), Eh = c(-0.5, 0.8, res), iprotein = ip)
  # Set up groups for affinity aggregation:
  # 3 Form I/(II)/III and 3 Form I proteins
  groups <- list(pre = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), post = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  amean <- agg.affinity(aout, groups = groups)
  diagram(amean, lwd = 2, col = 3, names = "", xlab = "pH", ylab = axis.label("Eh"), balance = 1)
  text(6, -0.18, "Higher mean affinity\nfor Form I/(II)/III proteins", srt = -21)
  text(6.5, 0.1, "Higher mean affinity\nfor Form I proteins", srt = -21)
  title("Groupwise Rubiscos", font.main = 1)
  label.figure("B", cex = 1.5, font = 2, yfrac = 0.92)

  # Panel C: Groupwise stability boundary (logfO2-pH diagram)
  par(mar = c(3, 3.5, 2.5, 3))
  plot_stability("rubisco_old", plot_names = FALSE, col = 3, O2lim = c(-70, -58))
  text(6.5, -64, "Higher mean affinity\nfor Form I/(II)/III proteins")
  text(5.5, -60, "Higher mean affinity\nfor Form I proteins")
  add_Eh7_axis()
  title("Groupwise Rubiscos", font.main = 1)
  label.figure("C", cex = 1.5, font = 2, yfrac = 0.92)

  if(pdf) dev.off()

}

# Figure S3: Relative stability of Rubisco stages (Kaçar et al., 2017 and Amritkar et al., 2025 individual datasets)  20260722
genoGOE_S3 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure_S3.pdf", width = 10, height = 4)
  par(mfrow = c(1, 2))
  # Kaçar et al., 2017
  plot_stability("rubisco_KHAB17_1_2", plot_names = FALSE, O2lim = c(-72.5, -40), col = 3)
  plot_stability("rubisco_KHAB17_2_3", plot_names = FALSE, add = TRUE, col = 3)
  text(7, -54, "Stage 1")
  text(7, -49, "Stage 2")
  text(6, -70, "Stage 2")
  text(6, -66, "Stage 3")
  abline(h = -62, lty = 2, col = 8)
  text(7, -61, "Upper limit of Fig. 6B", font = 3)
  label.figure("A", cex = 1.5, font = 2, yfrac = 0.94)
  title("Kaçar et al. (2017)", font.main = 1)
  # Amritkar et al., 2025
  plot_stability("rubisco_ACK25", O2lim = c(-72.5, -40), col = 3)
  abline(h = -62, lty = 2, col = 8)
  label.figure("B", cex = 1.5, font = 2, yfrac = 0.94)
  title("Amritkar et al. (2025)", font.main = 1)
  if(pdf) dev.off()
}

