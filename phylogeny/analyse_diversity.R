##### Plotting
library(here)
source(here('helper_functions.R'))

 

deduplicated_genus_tree = ape::read.tree(file.path('inputs','prepared_final_smbtree.tre'))
genus_abundance_diversity_data = read.csv(file.path('..','diversity_metrics','outputs','genus_level_pathway_diversity_information.csv'))
genus_distance_diversity_data = read.csv(file.path('..','diversity_metrics','outputs','genus_level_distance_diversity_information.csv'))
all_data = merge(genus_abundance_diversity_data,genus_distance_diversity_data, by= 'Genus')

significant_vars= c('FAD_minmax', 'G_minmax', 'H_minmax', 'Hbc_minmax', 'J_minmax', 'MFAD_minmax', 'APWD_minmax')
rename = significant_vars
data_to_use = all_data[c('Genus',significant_vars)]
labelled_tree = get_subset_of_tree_from_genera_in_data(data_to_use,deduplicated_genus_tree)
colnames(data_to_use) <- c('Genus',rename)
data_to_use = get_matching_genus_labels(labelled_tree,data_to_use)[c('label',rename)]
heatmap_plot(labelled_tree,data_to_use,'Diversity','diversity_plot.jpg',plotwidth = 6,plotheight = 5)
heatmap_plot(labelled_tree,data_to_use,'Diversity','diversity_plot_with_tips.jpg',plotwidth = 10,plotheight = 25, tipnames=TRUE)

calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'Hbc')
calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'H')
calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'J')
calculate_signal(deduplicated_genus_tree,genus_abundance_diversity_data,'G')

calculate_signal(deduplicated_genus_tree,genus_distance_diversity_data,'APWD')
calculate_signal(deduplicated_genus_tree,genus_distance_diversity_data,'FAD')
calculate_signal(deduplicated_genus_tree,genus_distance_diversity_data,'MFAD')


