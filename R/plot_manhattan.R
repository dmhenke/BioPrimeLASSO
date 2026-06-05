#' Manhattan plots of BioPrimeLASSO results
#'
#' Generates three Manhattan-style plots showing the correlation of each gene's
#' omic feature with the target dependency score, overlaid with non-zero
#' coefficients from the baseline and bio-primed LASSO models.
#'
#' @param gene A character string giving the gene symbol of interest (e.g., \code{"EGFR"}).
#' @param resIn Path to an \code{.RData} file containing a \code{results_omic} list
#'   as returned by \code{\link{bplasso}}, with an additional \code{cor2score} element.
#' @param subplotChr Integer chromosome number for the zoomed subplot. If \code{NA}
#'   (default), the chromosome of \code{gene} is used automatically.
#' @param dependency A data frame of dependency scores with cell lines as rows
#'   and gene symbols as column names.
#' @param gene_info A data frame of gene annotations with columns:
#'   \code{hgnc_symbol}, \code{chromosome_name}, and \code{start_position}.
#' @param dir_save Character string giving the directory path where the three
#'   PDF plots will be written.
#'
#' @return Called for side effects. Writes three PDF files to \code{dir_save}
#'   and renders the labeled full-genome plot to the active graphics device.
#' @export
#'
#' @examples
#' \dontrun{
#'   plot_manhattan(
#'     gene = "EGFR",
#'     resIn = "./EGFR_demeter2_CNV.RData",
#'     subplotChr = 11,
#'     dependency = demeter2,
#'     gene_info = gene_info,
#'     dir_save = "./Outputs/Graphics/")
#' }
plot_manhattan <- function(gene, resIn, subplotChr = NA, dependency, gene_info, dir_save){
  load(resIn)
  y <- dependency[, gene]
  y <- y[!is.na(y)]
  correl <- results_omic$cor2score
  aframe <- data.frame(
    gene = names(correl),
    gene_info[match(names(correl), gene_info$hgnc_symbol), ],
    correl
  )
  aframe <- aframe[order(aframe$chromosome_name, aframe$start_position), ]
  aframe$rank <- seq_len(nrow(aframe))

  subm <- results_omic
  subm_betas <- subm$betas[match(aframe$gene, rownames(subm$betas)), ]

  subm_betas$betas_pen[subm_betas$betas_pen == 0] <- NA
  subm_betas$betas[subm_betas$betas == 0] <- NA

  aframe$betas_pen <- subm_betas$betas_pen
  aframe$betas <- subm_betas$betas
  aframe$betalogic <- apply(cbind(aframe$betas_pen, aframe$betas), 1, function(x){
    if (is.na(x[1]) & is.na(x[2])) NA
    else if (!is.na(x[1]) & is.na(x[2])) "beta_pen"
    else if (!is.na(x[1]) & !is.na(x[2])) "both"
    else "beta"
  })
  aframe <- aframe[!is.na(aframe$chromosome_name), ]
  aframe$betasize <- aframe$betas
  aframe$betasize[!is.na(aframe$betas_pen)] <- aframe$betas_pen[!is.na(aframe$betas_pen)]
  aframe$betasize <- abs(aframe$betasize)

  lab_x <- "Gene ordered by genomic coordinate"
  lab_y <- paste0("Correlation (r)\n", gene, " Dependency with 'omic")

  g_full <- ggplot2::ggplot(aframe,
                   ggplot2::aes(rank, correl, color = chromosome_name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "red") +
    ggplot2::geom_point(size = 0.3) +
    ggplot2::scale_color_manual(
      guide = "none",
      breaks = c(1:22, "beta_pen", "both", "beta"),
      values = c(rep(c("black", "grey"), 11), "darkblue", "purple", "darkred")) +
    ggplot2::labs(x = lab_x, y = lab_y) +
    ggplot2::theme_classic()

  g_fullLab <- g_full +
    ggrepel::geom_text_repel(
      min.segment.length = 0,
      force = 1,
      direction = "both",
      max.overlaps = 100,
      max.time = 0.3,
      max.iter = 1e5,
      data = aframe[!is.na(aframe$betalogic), ],
      ggplot2::aes(label = gene, color = betalogic, size = betasize)) +
    ggplot2::scale_size(guide = "none")

  ggplot2::ggsave(
    filename = file.path(dir_save, paste0("PlotManhattan_", gene, ".pdf")),
    plot = g_full, width = 10, height = 5)
  ggplot2::ggsave(
    filename = file.path(dir_save, paste0("PlotManhattan_", gene, "_labs.pdf")),
    plot = g_fullLab, width = 10, height = 5)

  # Subplot of target chromosome
  if (is.na(subplotChr)) {
    which_chrome <- aframe[aframe$gene == gene, "chromosome_name"]
  } else {
    which_chrome <- subplotChr
  }
  which_chromePos <- which(aframe$chromosome_name == which_chrome)
  which_chrome_minX <- min(which_chromePos)
  which_chrome_maxX <- max(which_chromePos)
  which_chrome_minY <- min(aframe[aframe$chromosome_name == which_chrome, "correl"])
  which_chrome_maxY <- max(aframe[aframe$chromosome_name == which_chrome, "correl"])

  g_chrGene <- g_fullLab +
    ggplot2::coord_cartesian(
      xlim = c(which_chrome_minX, which_chrome_maxX),
      ylim = c(which_chrome_minY, which_chrome_maxY)) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = paste0("Chromosome ", which_chrome))

  ggplot2::ggsave(
    filename = file.path(dir_save, paste0("PlotManhattan_", gene, "_chr", which_chrome, ".pdf")),
    plot = g_chrGene, width = 4, height = 2)

  plot(g_fullLab)
}
