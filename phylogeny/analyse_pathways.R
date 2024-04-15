
library(here)
source(here('helper_functions.R'))
np_pathways = c('Terpenoids', 'Fatty_acids', 'Polyketides', 'Carbohydrates', 'Amino_acids_and_Peptides', 'Shikimates_and_Phenylpropanoids',
               'Alkaloids')
deduplicated_genus_tree = ape::read.tree(file.path('inputs','prepared_final_smbtree.tre'))
genus_level_data = read.csv(file.path('..','collect_compound_data','outputs','genus_level_pathway_data.csv'))
labelled_tree = get_subset_of_tree_from_genera_in_data(genus_level_data,deduplicated_genus_tree)


for (pathway in np_pathways){
  
  calculate_signal(deduplicated_genus_tree,genus_level_data,paste('mean_identified_as',pathway,sep='_'))
  calculate_signal(deduplicated_genus_tree,genus_level_data,paste('norm_mean_identified_as',pathway,sep='_'))
}


# Vars that are siginficant in both standard and normed cases
significant_vars= c('mean_identified_as_Polyketides',  'mean_identified_as_Alkaloids','mean_identified_as_Terpenoids','mean_identified_as_Shikimates_and_Phenylpropanoids')
rename = c("Polyketides", "Alkaloids", "Terpenoids", 'Shikimates_and_Phenylpropanoids')
data_to_use = genus_level_data[c('Genus',significant_vars)]
colnames(data_to_use) <- c('Genus',rename)
data_to_use = get_matching_genus_labels(labelled_tree,data_to_use)[c('label',rename)]
heatmap_plot(labelled_tree,data_to_use,"Proportion Identified\nAs Pathway",'pathway_plot.jpg',plotwidth = 6,plotheight = 5)
heatmap_plot(labelled_tree,data_to_use,"Proportion Identified\nAs Pathway",'pathway_plot_with_tips.jpg',plotwidth = 10,plotheight = 25, tipnames=TRUE)
