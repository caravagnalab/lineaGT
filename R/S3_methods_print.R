#' Print method
#'
#' @description add
#'
#' @param x An object of class mvnmm
#' @param ... Default extra paramaters
#'
#' @return Prints to screen information regarding the fitted object.
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

  tp = x %>% get_tp_to_int() %>% unlist() %>% sort() %>% names()
  if (purrr::is_empty(tp)) tp = x %>% get_timepoints()

  cat("\n")
  cli::cli_text(clisymbols::symbol$arrow_right, " Lineages: {.field {get_lineages(x) %>% sort()}}.")
  cli::cli_text(clisymbols::symbol$arrow_right, " Timepoints: {.field {tp}}.")
  cli::cli_text(clisymbols::symbol$arrow_right, " Number of Insertion Sites: {.field {x$data.shape[1]}}.")

  pi = x %>% get_weights() %>% round(2) %>% sort(decreasing = TRUE)
  n_IS = x %>% get_ISs() %>% sort(decreasing=T)

  cli::cli_h3("Optimal IS model with {.field k = {pi %>% length}}.")
  cat("\n")

  for(cluster in names(n_IS))
  {
    starting =
      sprintf("%25s", paste0(crayon::yellow(cluster),
                             " (",
                             # pi[cluster] * 100,
                             # "% - ",
                             n_IS[cluster],
                             " ISs)")
              )

    # Gaussian means
    mus = (x %>% get_mean())[cluster, ] %>% round(0)

    n_chars = x %>% get_mean() %>% round(0) %>% as.character() %>% nchar() %>% max()
    mus = sapply(mus, function(x) sprintf(paste0("%", n_chars , 's'), x))

    lins = get_lineages(x) %>% sort()
    # tmp = x %>% get_tp_to_int() %>% unlist() %>% sort() %>% names()
    tmp = tp

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
