plot_losses = function(x, train=FALSE) {
  losses = get_losses(x, x$runs, train=train)

  if (!train)
    return(losses %>%
             ggplot() +
             geom_line(aes(x=index, y=losses)) +
             scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
             ylab("ELBO") + xlab("Iterations") +
             my_ggplot_theme()
    )

  losses = losses %>%
    dplyr::mutate(losses=purrr::map(losses, ~data.frame(losses=.x, index=seq_along(.x)))) %>%
    tidyr::unnest(losses) %>%
    dplyr::mutate(K=factor(K, levels=K %>% unique()), run=factor(run, levels=run %>% unique())) %>%
    tidyr::separate("id", into=c("id.K", "id.run"), sep="[.]")

  return(losses %>%
           ggplot() +
           # geom_point(aes(x=index, y=losses, color=K), size=.5) +
           geom_line(aes(x=index, y=losses, color=K, group=id.K), size=.5) +
           facet_wrap(~run) +
           scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
           ylab("ELBO") + xlab("Iterations") +
           my_ggplot_theme()
         )

}

plot_gradient_norms = function(x) {
  grads = x$runs %>% get_gradient_norms() %>%
    mutate(grad_norm=purrr::map(grad_norm, ~data.frame(grad_norm=.x, index=seq_along(.x)))) %>%
    tidyr::unnest(grad_norm) %>%
    mutate(K=factor(K, levels=K %>% unique()), run=factor(run, levels=run %>% unique())) %>%
    tidyr::separate("id", into=c("id.K", "id.run"), sep="[.]")

  return(grads %>%
           ggplot() +
           geom_line(aes(x=index, y=grad_norm, color=K, group=id.K), size=.5) +
           facet_grid(rows=vars(param), cols=vars(run)) +
           scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
           ylab("Gradient Norms") + xlab("Iterations") +
           my_ggplot_theme()
         )
}


plot_IC = function(x) {
  ic = x$runs$IC %>%
    tidyr::separate("id", into=c("id.K", "id.run"), sep="[.]")

  return(ic %>%
           ggplot() +
           # geom_point(aes(x=factor(K), y=mean_val, color=as.factor(run))) +
           # geom_errorbar(aes(x=factor(K), ymin=min_val, ymax=max_val)) +
           geom_point(aes(x=factor(K), y=value, color=as.factor(run))) +
           geom_line(aes(x=factor(K), y=value, color=as.factor(run), group=run)) +
           facet_wrap(~method, scales="free_y") +
           xlab("K") + ylab("Value") + labs(color="") +
           scale_x_discrete(limits=factor(ic$K %>% unique %>% sort)) +
           my_ggplot_theme()
         )
}


