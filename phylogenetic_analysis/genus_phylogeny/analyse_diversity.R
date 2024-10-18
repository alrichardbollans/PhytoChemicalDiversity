##### Plotting
library(here)
source(here('helper_functions.R'))

 

deduplicated_genus_tree = ape::read.tree(genus_tree_path)
all_data = read.csv(file.path('..','..','collect_and_compile_data','get_diversity_metrics', 'outputs', 'group_data', 'Genus_transformed.csv'))


vars_to_plot = c( 'FAD', 'MFAD', 'APWD', 'H', 'Hbc','G')
data_to_use = all_data[c('Assigned_group','N',vars_to_plot)]
colnames(data_to_use) <- c('Genus','N',vars_to_plot)
data_to_use <- na.omit(data_to_use)
labelled_tree = get_subset_of_tree_from_genera_in_data(data_to_use,deduplicated_genus_tree)

data_to_plot = get_matching_genus_labels(labelled_tree,data_to_use)[c('label',vars_to_plot)]
heatmap_plot(labelled_tree,data_to_plot,'Diversity','diversity_plot.jpg',plotwidth = 5,plotheight = 4)
heatmap_plot(labelled_tree,data_to_plot,'Diversity','diversity_plot_with_tips.jpg',labelsize=1.5,plotwidth = 10,plotheight = 10, tipnames=TRUE)

calculate_signal(deduplicated_genus_tree,data_to_use,'Hbc')
calculate_signal(deduplicated_genus_tree,data_to_use,'H')
calculate_signal(deduplicated_genus_tree,data_to_use,'G')

calculate_signal(deduplicated_genus_tree,data_to_use,'APWD')
calculate_signal(deduplicated_genus_tree,data_to_use,'FAD')
calculate_signal(deduplicated_genus_tree,data_to_use,'MFAD')

calculate_signal(deduplicated_genus_tree,data_to_use,'N')

