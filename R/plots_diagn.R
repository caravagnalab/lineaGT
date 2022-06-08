plot_losses = function(x) {
  losses = x %>% get_losses() %>%
    mutate(losses=purrr::map(losses, ~data.frame(losses=.x, index=seq_along(.x)))) %>%
    tidyr::unnest(losses) %>%
    mutate(K=factor(K, levels=K %>% unique()), run=factor(run, levels=run %>% unique()))

  return(losses %>%
           ggplot() +
           geom_line(aes(x=index, y=losses, color=K)) +
           facet_wrap(~run) +
           scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
           ylab("Losses") + xlab("Iterations")
         )

}

plot_gradient_norms = function(x) {
  grads = x %>% get_gradient_norms() %>%
    mutate(grad_norm=purrr::map(grad_norm, ~data.frame(grad_norm=.x, index=seq_along(.x)))) %>%
    tidyr::unnest(grad_norm) %>%
    mutate(K=factor(K, levels=K %>% unique()), run=factor(run, levels=run %>% unique()))

  return(grads %>%
           ggplot() +
           geom_line(aes(x=index, y=grad_norm, color=K)) +
           facet_wrap(param~run) +
           scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
           ylab("Gradient Norms") + xlab("Iterations")
         )
}

plot_IC = function(x) {
  ic = x %>% get_IC()

  return(ic %>%
           ggplot() +
           geom_line(aes(x=K, y=value, color=method)) +
           facet_wrap(method~run, scales="free_y", nrow=ic$method %>% unique() %>% length()) +
           xlab("K") + ylab("Value")
         )

}


