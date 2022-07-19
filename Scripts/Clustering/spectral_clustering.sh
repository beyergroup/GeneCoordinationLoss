#!/bin/bash

# Laplacian, |weight|
Rscript --vanilla spectral_clustering.R "stabsel" "network_largest_cc" "abs" "Laplacian" >> ../../Logs/Laplacian_absweights_spectral_clustering.out

# Laplacian, no weights
Rscript --vanilla spectral_clustering.R "stabsel" "network_largest_cc" "none" "Laplacian" >> ../../Logs/Laplacian_noweights_spectral_clustering.out

# Adjacency, |weight|
Rscript --vanilla spectral_clustering.R "stabsel" "network_largest_cc" "abs" "Adjacency" >> ../../Logs/Adjacency_absweights_spectral_clustering.out

# Adjacency, no weights
Rscript --vanilla spectral_clustering.R "stabsel" "network_largest_cc" "none" "Adjacency" >> ../../Logs/Adjacency_noweights_spectral_clustering.out
