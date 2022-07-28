plot_losses = function(x, train=FALSE) {
  losses = x %>% get_losses(train=train)

  if (!train)
    return(losses %>%
             ggplot() +
             geom_line(aes(x=index, y=losses)) +
             # scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
             ylab("ELBO") + xlab("Iterations")
    )

  losses = losses %>%
    mutate(losses=purrr::map(losses, ~data.frame(losses=.x, index=seq_along(.x)))) %>%
    tidyr::unnest(losses) %>%
    mutate(K=factor(K, levels=K %>% unique()), run=factor(run, levels=run %>% unique()))

  return(losses %>%
           ggplot() +
           geom_point(aes(x=index, y=losses, color=K)) +
           facet_wrap(~run) +
           scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
           ylab("ELBO") + xlab("Iterations") +
           my_ggplot_theme()
         )

}

plot_gradient_norms = function(x) {
  grads = x %>% get_gradient_norms() %>%
    mutate(grad_norm=purrr::map(grad_norm, ~data.frame(grad_norm=.x, index=seq_along(.x)))) %>%
    tidyr::unnest(grad_norm) %>%
    mutate(K=factor(K, levels=K %>% unique()), run=factor(run, levels=run %>% unique()))

  return(grads %>%
           ggplot() +
           geom_point(aes(x=index, y=grad_norm, color=K), size=.5) +
           facet_grid(rows=vars(param), cols=vars(run)) +
           scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
           ylab("Gradient Norms") + xlab("Iterations")
         )
}

plot_IC = function(x) {
  ic = x %>% get_IC()

  # ic %>%
  #   ggplot() +
  #   geom_point(aes(x=K, y=value, color=method)) +
  #   facet_grid(method~run, scales="free_y") +
  #   xlab("K") + ylab("Value") +
  #   scale_x_discrete(limits = ic$K %>% unique %>% sort)
  #
  # return(ic %>%
  #          ggplot() +
  #          geom_point(aes(x=K, y=value, color=method)) +
  #          facet_grid(method~run, scales="free_y", nrow=ic$method %>% unique() %>% length()) +
  #          xlab("K") + ylab("Value")
  #        )

  ic %>%
    ggplot() +
    geom_line(aes(x=K, y=value, color=method)) +
    geom_point(aes(x=K, y=value, color=method)) +
    facet_wrap(~method, scales="free_y") +
    xlab("K") + ylab("Value") +
    scale_x_discrete(limits = ic$K %>% unique %>% sort) +
    my_ggplot_theme()
}


