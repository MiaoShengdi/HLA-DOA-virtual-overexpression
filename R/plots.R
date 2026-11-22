## ============================================================================
##  R/plots.R
##
##  Figures for the differential-regulation output. Aesthetics are kept plain
##  (white background, two-colour coding) so the figures read well in print.
## ============================================================================

.sig_colour  <- "#B23A2E"
.ns_colour   <- "#AEBABB"
.bar_colour  <- "#2F5D62"

#' Top genes by network perturbation score
plot_top_genes <- function(diff, top_n = 20) {

  d <- diff[order(-diff$FC), , drop = FALSE]
  d <- d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
  d$gene <- factor(d$gene, levels = d$gene[order(d$FC)])

  ggplot2::ggplot(d, ggplot2::aes(x = gene, y = FC)) +
    ggplot2::geom_col(fill = .bar_colour, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Gene", y = "Network perturbation score (FC)") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
}

#' Displacement against significance
plot_volcano <- function(diff, fdr = 0.05, n_label = 10) {

  d <- diff
  ## p-values can underflow to zero for a strongly displaced gene; floor them
  ## so that -log10() stays finite, and clip the axis at the 99.9% quantile
  d$neg_log10_fdr <- -log10(pmax(d$p.adj, 1e-300))
  y_max <- stats::quantile(d$neg_log10_fdr, 0.999, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  d$status <- ifelse(d$p.adj < fdr, "Significant", "Not significant")

  lab <- d[d$p.adj < fdr & d$neg_log10_fdr <= y_max, , drop = FALSE]
  if (nrow(lab) > 0) {
    lab <- lab[order(-lab$FC), , drop = FALSE]
    lab <- lab[seq_len(min(n_label, nrow(lab))), , drop = FALSE]
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = Z, y = neg_log10_fdr,
                                       colour = status)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = c("Significant" = .sig_colour,
                                            "Not significant" = .ns_colour)) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed",
                        linewidth = 0.4, colour = .sig_colour) +
    ggplot2::coord_cartesian(ylim = c(0, y_max)) +
    ggplot2::labs(x = "Z (network displacement)",
                  y = expression(-log[10] * "(FDR)"), colour = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "none")

  if (nrow(lab) > 0) {
    p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene),
                                      size = 2.8, max.overlaps = 30,
                                      colour = "#3A3A3A", seed = 1)
  }
  p
}

#' Proportion of significantly displaced genes
plot_proportion <- function(diff, fdr = 0.05) {

  status <- factor(ifelse(diff$p.adj < fdr, "Significant", "Not significant"),
                   levels = c("Significant", "Not significant"))
  tab <- as.data.frame(table(status), stringsAsFactors = FALSE)
  names(tab) <- c("status", "n")
  tab$label <- paste0(round(100 * tab$n / sum(tab$n), 1), "%")

  ggplot2::ggplot(tab, ggplot2::aes(x = "", y = n, fill = status)) +
    ggplot2::geom_col(width = 1, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       position = ggplot2::position_stack(vjust = 0.5),
                       size = 3.2, colour = "#2B2B2B") +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::scale_fill_manual(values = c("Significant" = .sig_colour,
                                          "Not significant" = "#D8DEDE")) +
    ggplot2::labs(fill = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right")
}

#' Write a figure as PNG (300 dpi) and PDF
save_figure <- function(plot, path_stem, width = 6, height = 5) {
  ggplot2::ggsave(paste0(path_stem, ".png"), plot, width = width,
                  height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(paste0(path_stem, ".pdf"), plot, width = width,
                  height = height, bg = "white")
}
