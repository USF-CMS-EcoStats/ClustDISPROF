#' Plot dendrogram of a DISPROF-based custer analysis
#'
#' This function plots the dendrogram based on the result from f_disprof_clust
#' analysis.
#'
#' @param result Output list from the 'f_disprof_clust' function
#' @param orient Plot orientation
#' \itemize{
#'  \item Vertical ('V')
#'  \item Horizontal ('H')
#' }
#'
#' @import grDevices
#'
#' @importFrom magrittr %>%
#' @importFrom dendextend set color_branches
#' @importFrom colorspace darken
#'
#' @export

f_disprof_clustPlot <- function(result,
                                orient = 'V'){

dnd       <- as.dendrogram(result$Z)
colors    <- darken(rainbow(n = length(unique(result$grp))), 0.2)
temp_cols <- colors[result$grp]
temp_cols <- temp_cols[order.dendrogram(dnd)]
temp_cols <- factor(temp_cols, unique(temp_cols))

dnd <- dnd %>%
  color_branches(clusters = as.numeric(temp_cols),
                 col     = levels(temp_cols)) %>%
  set("labels_colors", as.character(temp_cols)) %>%
  # Adjust line-weight of branches
  set("branches_lwd", 1.5) %>%
  # Adjust size of labels
  set("labels_cex", 0.75)

  if (orient == 'V') {
    plot(dnd,
         horiz = FALSE,
         main = "DISPROF-based Cluster Analysis")
  } else if (orient == "H") {
    plot(dnd,
         horiz = TRUE,
         main = "DISPROF-based Cluster Analysis")
  }
}
