## Predictability analysis

# First, predict expression values using the networks in GTEx subsets:
predict_expression.R

# Then, compute correlation between predicted and observed values:
compute_network_correlations.R

# This is done per tissue and for the cross-tissue (all subsets of the same size).

# Then, analyse overlap of poorly predicted genes across tissues:
compute_plot_poorly_predicted.R

# Characterize the poorly predicted genes based on dataset- and network-specific features:
characterize_poorly_predicted.R

# Compare predictability with and without partial correlations:
plot_predictability_comparison_GTEx_subsets.R
plot_predictability_stabsel_on_pcclasso.R



## Topology analysis

# Look at in- vs out-going edges in stabsel_pcclasso network:
in_out_edges.R

# Remove unbalanced edges (should be directed but are not):
remove_unbalanced_edgepairs.R

# Plot centrality metrics:
topology.R



## Predictability changes with age

# Compute difference between predicted by network and observed, for each gene in each sample:
compute_network_delta.R

# Limma for both expression and predictability changes with age:
compute_age_limma_GTEx.R

# Plot expression vs predictability:
plot_age_cordeltas_logFCs_tissues.R

# GO enrichment in genes with significant predictabilty changes with age:
expression_predictability_FCs.R