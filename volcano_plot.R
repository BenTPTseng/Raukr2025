#' Volcano plot
#'
#' @param x a data frame containing CpG location, strand information and statistics
#' @param qvalue_cutoff 0.05
#' @param meth.diff_cutoff 25
#' @param hypermeth_col "red3"
#' @param hypometh_col "royalblue3"
#' @param common_col "seashell3"
#'
#' @return ggplot object
#' @export
#'
volcano_plot <- function(x,
                         qvalue_cutoff    = 0.05,
                         meth.diff_cutoff = 25,
                         hypermeth_col    = "red3",
                         hypometh_col     = "royalblue3",
                         common_col       = "seashell3"){

    tmp <- x

    hypermeth <- ((tmp$qvalue < qvalue_cutoff) & (tmp$meth.diff >  meth.diff_cutoff))
    hypometh  <- ((tmp$qvalue < qvalue_cutoff) & (tmp$meth.diff < -meth.diff_cutoff))

    tmp$Group <- rep("common", times = nrow(tmp))

    tmp[hypermeth, "Group"]  <- "hypermeth"
    tmp[hypometh,  "Group"]  <- "hypometh"

    tmp$Group <- factor(tmp$Group, levels = c("hypermeth", "common", "hypometh"))

    ggplot(data = tmp,
           aes(x = meth.diff, y = -log10(qvalue), colour = Group)) +
      xlab(label = "Methylation percentage difference") +
      ylab(label = "-log10(Adjusted p-value)") +
      geom_point(alpha = 0.8) +
      xlim(c(-100, 100)) +
      ylim(c(-0, 100)) +
      geom_hline(yintercept = -log10(qvalue_cutoff), linetype = "dashed") +
      geom_vline(xintercept = c(meth.diff_cutoff, -meth.diff_cutoff), linetype = "dashed") +
      annotate("text", x =  75 , y = 90, label = paste0("n = ", sum(hypermeth))) +
      annotate("text", x = -75 , y = 90, label = paste0("n = ", sum(hypometh))) +
      annotate("text", x = max(tmp$meth.diff) + 10, y = -log10(qvalue_cutoff) - 1, label = "p=0.05") +
      theme_classic() +
      scale_color_manual(values = c(hypermeth_col, common_col, hypometh_col),
                         labels = c("Hypermethylated CpGs", "Non-significant CpGs", "Hypomethylated CpGs")) +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
  }
