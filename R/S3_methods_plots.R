#' Mullerplot
#'
#' @description Mullerplot showing the longitudinal clonal evolution per lineage.
#'
#' @param x An object of class mvnmm
#' @param ... Default extra parameters
#'
#' @return A ggplot object with the mullerplot for the fitted object.
#'
#' @exportS3Method plot mvnmm
#' @export plot.mvnmm


plot.mvnmm = function(x, ...) {
  plot_mullerplot(x, ...)
}

my_ggplot_theme = function(cex=1, legend.pos="right", text_size=9) {
  theme_light(base_size=10 * cex) +
    theme(legend.position=legend.pos, legend.key.size=unit(.3*cex, "cm"),
          panel.background=element_rect(fill='white'), legend.key.width=unit(.5*cex, 'cm'),
          text=element_text(size=text_size))
}
