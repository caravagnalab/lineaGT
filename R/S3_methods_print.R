#' Print method
#'
#' @description add
#'
#' @param x An object of class mvnmm
#' @param ... Default extra paramaters
#'
#' @return
#'
#' @exportS3Method print mvnmm
#' @export print.mvnmm

print.mvnmm = function(x, ...) {

  stopifnot(inherits(x, "mvnmm"))

  pylog = reticulate::py_discover_config()$python

  cli::cli_rule(
    left = paste(
      crayon::bgYellow(crayon::black("[ lineaGT ] "))),
    right = paste(crayon::bold("Python:"), pylog %>% crayon::silver())
    )

  which_conda_env()
  cat("\n")
  cli::cli_text(clisymbols::symbol$arrow_right, " Lineages: {.field {get_lineages(x) %>% sort()}}.")
  cli::cli_text(clisymbols::symbol$arrow_right, " Timepoints: {.field {get_timepoints(x)}}.")

  pi = x %>% lineaGT::get_weights() %>% round(2) %>% sort(decreasing = TRUE)

  cli::cli_h3("Optimal IS model with {.field k = {pi %>% length}}.")
  cat("\n")

  for(cluster in names(pi))
  {
    starting =
      sprintf("%25s", paste0(
        crayon::yellow(cluster),
        ' (', pi[cluster] * 100, '%)'
      ))

    # Gaussian means
    mus = (x %>% get_mean())[cluster, ] %>% round(0)

    n_chars = x %>% get_mean() %>% round(0) %>% as.character() %>% nchar() %>% max()
    mus = sapply(mus, function(x) sprintf(paste0("%", n_chars , 's'), x))

    lins = get_lineages(x) %>% sort()
    tmp = get_timepoints(x)

    inliners = paste0(
      crayon::blue(lins),
      ' [',
        lapply(lins, function(l) paste('cov', tmp, l, sep = '.')) %>%
          sapply(function(w) paste0(mus[w], collapse = ', ')),
      ']'
    )

    paste(starting, ':', inliners %>% paste(collapse = '; '), "\n") %>%
      cat
  }


  return(x)
}
