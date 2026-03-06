# Package index

## Fit function

Main function to perform the fit of the model.

- [`fit()`](caravagnalab.github.io/lineaGT/reference/fit.md) :

  Creates an object of class `mvnmm`.

- [`fit_mutations()`](caravagnalab.github.io/lineaGT/reference/fit_mutations.md)
  : Fit the mutations clustering

- [`fit_phylogenies()`](caravagnalab.github.io/lineaGT/reference/fit_phylogenies.md)
  : Fit the phylogenetic trees

- [`fit_growth_rates()`](caravagnalab.github.io/lineaGT/reference/fit_growth_rates.md)
  : Infer growth rates for each clone and subclone.

- [`filter_dataset()`](caravagnalab.github.io/lineaGT/reference/filter_dataset.md)
  : Filters the input dataset.

## Getter functions

Functions to extract the main elements of the fitted object.

- [`get_lineages()`](caravagnalab.github.io/lineaGT/reference/get_lineages.md)
  : Extract the data lineages.
- [`get_dimensions()`](caravagnalab.github.io/lineaGT/reference/get_dimensions.md)
  : Extract the model dimensions.
- [`get_timepoints()`](caravagnalab.github.io/lineaGT/reference/get_timepoints.md)
  : Extract the data timepoints.
- [`get_cov_dataframe()`](caravagnalab.github.io/lineaGT/reference/get_cov_dataframe.md)
  : Retrieve the coverage dataframe.
- [`get_vaf_dataframe()`](caravagnalab.github.io/lineaGT/reference/get_vaf_dataframe.md)
  : Retrieve the mutations dataframe.
- [`get_labels()`](caravagnalab.github.io/lineaGT/reference/get_labels.md)
  : Extract the observations labels.
- [`get_unique_labels()`](caravagnalab.github.io/lineaGT/reference/get_unique_labels.md)
  : Extract the list of unique observations labels.
- [`get_unique_muts_labels()`](caravagnalab.github.io/lineaGT/reference/get_unique_muts_labels.md)
  : Retrieve the list of unique labels of mutation clusters.
- [`get_mean()`](caravagnalab.github.io/lineaGT/reference/get_mean.md) :
  Extract the estimated mean parameters.
- [`get_sigma()`](caravagnalab.github.io/lineaGT/reference/get_sigma.md)
  : Extract the estimated variance parameters.
- [`get_covariance_Sigma()`](caravagnalab.github.io/lineaGT/reference/get_covariance_Sigma.md)
  : Extract the estimated covariance matrices.
- [`get_covariance_Cholesky()`](caravagnalab.github.io/lineaGT/reference/get_covariance_Cholesky.md)
  : Extract the estimated Cholesky matrices, used to factorise the
  covariance matrix.
- [`get_weights()`](caravagnalab.github.io/lineaGT/reference/get_weights.md)
  : Extract the estimated mixing proportions.
- [`get_z_probs()`](caravagnalab.github.io/lineaGT/reference/get_z_probs.md)
  : Extract the estimated posterior probabilities.
- [`get_ISs()`](caravagnalab.github.io/lineaGT/reference/get_ISs.md) :
  Get the number of ISs per cluster.
- [`estimate_n_pops()`](caravagnalab.github.io/lineaGT/reference/estimate_n_pops.md)
  : Function implemented to estimate the real number of clones in each
  cluster.

## Visualization functions

Functions to visualize the fit results.

- [`plot_scatter_density()`](caravagnalab.github.io/lineaGT/reference/plot_scatter_density.md)
  : 2D scatterplot and density
- [`plot_mixture_weights()`](caravagnalab.github.io/lineaGT/reference/plot_mixture_weights.md)
  : Barplot of the per-cluster mixture weights and number of ISs.
- [`plot_marginal()`](caravagnalab.github.io/lineaGT/reference/plot_marginal.md)
  : Histogram of the marginal distribution of each dimension
- [`plot_mullerplot()`](caravagnalab.github.io/lineaGT/reference/plot_mullerplot.md)
  : Muller plot
- [`plot_growth_regression()`](caravagnalab.github.io/lineaGT/reference/plot_growth_regression.md)
  : Visualize the regression given the infered growth rates.
- [`plot_growth_rates()`](caravagnalab.github.io/lineaGT/reference/plot_growth_rates.md)
  : Visualize the infered growth rates.
- [`plot_vaf()`](caravagnalab.github.io/lineaGT/reference/plot_vaf.md) :
  VAF 2D scatterplot
- [`plot_vaf_time()`](caravagnalab.github.io/lineaGT/reference/plot_vaf_time.md)
  : VAF over time
- [`plot_phylogeny()`](caravagnalab.github.io/lineaGT/reference/plot_phylogeny.md)
  : Clonal evolution trees
- [`plot_differentiation_tree()`](caravagnalab.github.io/lineaGT/reference/plot_differentiation_tree.md)
  : Visualize the number of subclones on the differentiation tree

## Training

Functions to visualize losses, Information Criteria and gradient norms
computed during training.

- [`plot_losses()`](caravagnalab.github.io/lineaGT/reference/plot_losses.md)
  : Function to plot the training losses.
- [`plot_IC()`](caravagnalab.github.io/lineaGT/reference/plot_IC.md) :
  Function to plot the Information Criteria computed during model
  selection.
- [`plot_gradient_norms()`](caravagnalab.github.io/lineaGT/reference/plot_gradient_norms.md)
  : Function to plot the gradients norms.

## S3 object methods

Print and plot methods for the fitted object.

- [`print(`*`<mvnmm>`*`)`](caravagnalab.github.io/lineaGT/reference/print.mvnmm.md)
  : Print method
- [`plot(`*`<mvnmm>`*`)`](caravagnalab.github.io/lineaGT/reference/plot.mvnmm.md)
  : Mullerplot

## Data

Example of the input datasets and of the fitted object.

- [`cov.df.example`](caravagnalab.github.io/lineaGT/reference/cov.df.example.md)
  : Example coverage data
- [`vaf.df.example`](caravagnalab.github.io/lineaGT/reference/vaf.df.example.md)
  : Example mutation data
- [`vaf.df.example`](caravagnalab.github.io/lineaGT/reference/x.example.md)
  : Example mutation data
