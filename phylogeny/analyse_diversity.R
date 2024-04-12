##### Plotting
library(here)
source(here('helper_functions.R'))

 

deduplicated_genus_tree = ape::read.tree(file.path('inputs','prepared_final_smbtree.tre'))
genus_abundance_diversity_data = read.csv(file.path('..','diversity_metrics','outputs','genus_level_pathway_diversity_information.csv'))
genus_distance_diversity_data = read.csv(file.path('..','diversity_metrics','outputs','genus_level_distance_diversity_information.csv'))

# TODO: Add species_richness although have different matching genera
# Need to do as separate heatmaps as the scale is different
significant_vars= c('N')
rename = c('N')
data_to_use = genus_distance_diversity_data[c('Genus',significant_vars)]
labelled_tree = get_subset_of_tree_from_genera_in_data(data_to_use,deduplicated_genus_tree)
colnames(data_to_use) <- c('Genus',rename)
data_to_use = get_matching_genus_labels(labelled_tree,data_to_use)[c('label',rename)]
heatmap_plot(labelled_tree,data_to_use,'Diversity','diversity_plot.jpg')

diversity_plot(deduplicated_genus_tree,genus_abundance_diversity_data,'bc_shannon')
diversity_plot(deduplicated_genus_tree,genus_abundance_diversity_data,'shannon')
diversity_plot(deduplicated_genus_tree,genus_abundance_diversity_data,'pielou')
diversity_plot(deduplicated_genus_tree,genus_abundance_diversity_data,'simpson')
diversity_plot(deduplicated_genus_tree,genus_distance_diversity_data,'N')
diversity_plot(deduplicated_genus_tree,genus_distance_diversity_data,'APWD')
diversity_plot(deduplicated_genus_tree,genus_distance_diversity_data,'FAD')
diversity_plot(deduplicated_genus_tree,genus_distance_diversity_data,'MFAD')

calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'bc_shannon')
calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'shannon')
calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'pielou')
calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'simpson')

calculate_signal(deduplicated_genus_tree,genus_distance_diversity_data,'APWD')
calculate_signal(deduplicated_genus_tree,genus_distance_diversity_data,'FAD')
calculate_signal(deduplicated_genus_tree,genus_distance_diversity_data,'MFAD')


