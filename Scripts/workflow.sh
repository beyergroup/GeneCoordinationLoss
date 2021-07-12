## Predictability analysis

# First, predict expression values using the networks in GTEx subsets:
predict_expression.R

# Then, compute correlation between predicted and observed values:
compute_network_correlations.R

This is done per tissue and for the cross-tissue (all subsets of the same size).

# Then, analyse overlap of poorly predicted genes across tissues:
compute_plot_poorly_predicted.R

# Characterize the poorly predicted genes based on dataset- and network-specific features:
characterize_poorly_predicted.R

# Compare predictability with and without partial correlations:
plot_predictability_comparison_GTEx_subsets.R
