x = readRDS("/home/elenab/dati_elenab/mutation_dynamics/lineaGT_data/murine_tp/HOMO_PGK/E2P7/HOMO_PGK.E2P7.fit.Rds")

x %>% estimate_mean_ISs()
x %>% estimate_n_pops()
x %>% get_ISs()

mp1 = x %>% plot_mullerplot(estimate_npops=F)
mp2 = x %>% plot_mullerplot(estimate_npops=T)
mix_w = x %>% plot_mixture_weights()

patchwork::wrap_plots(mp1, mp2, mix_w, design="AAACCC
                                               BBB###")
