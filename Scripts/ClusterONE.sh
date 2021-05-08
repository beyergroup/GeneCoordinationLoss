#!/bin/bash

# Run ClusterONE on Human network
cd ../Outputs/Human_Network/Topology/ClusterONE
for D in 0.3 0.5 0.7
do
  mkdir $D
  # java -jar ../../../../Resources/cluster_one-1.0.jar -F csv -d $D ../../undirected_edge_list.txt >> ${D}/clusters.csv
  java -jar ../../../../Resources/cluster_one-1.0.jar -F csv -d $D ../../undirected_abs_scale_weights_edge_list.txt >> ${D}/clusters_weighted.csv
done
cd ../../../../Scripts
