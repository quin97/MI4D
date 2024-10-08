ggrare <- function (physeq_object, step = 10, label = NULL, color = NULL, 
                    plot = TRUE, title = "default", parallel = FALSE, se = TRUE) 
{
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) {
    x <- t(x)
  }
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), 
        sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, , drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    }
    else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  }
  else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), 
                       "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = "Size", 
                                                        y = ".S", group = "Sample", color = color))
  p <- p + ggplot2::labs(x = "Number of Reads", 
                         y = "ASV Richness",
                         title = title)
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels, ggplot2::aes_string(x = "x", 
                                                                   y = "y", label = label, color = color), size = 5, show.legend = FALSE,
                                hjust = 3)
  }
  p <- p + ggplot2::geom_line(size = 1)
  if (se) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se", 
                                                      ymax = ".S + .se", color = NULL, fill = color), 
                                  alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
} 

