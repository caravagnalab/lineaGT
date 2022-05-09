#' 2D scatter and density plot
#'
#' @description 2D scatter and density plot for each pair of dimensions
#'
#' @param x An object of class mvnmm
#' @param ... Default extra paramaters
#'
#' @return
#'
#' @exportS3Method plot mvnmm
#' @export plot.mvnmm


plot.mvnmm = function(x, ...) {

}

my_ggplot_theme = function(cex=1, legend.pos="right") {
  theme_light(base_size=10 * cex) +
    theme(legend.position=legend.pos, legend.key.size=unit(.3*cex, "cm"),
          panel.background=element_rect(fill='white'), legend.key.width=unit(.5*cex, 'cm'),
          text=element_text(size=9))
}