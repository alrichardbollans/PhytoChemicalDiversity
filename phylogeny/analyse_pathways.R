
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggnewscale)

np_pathways = c('Terpenoids', 'Fatty_acids', 'Polyketides', 'Carbohydrates', 'Amino_acids_and_Peptides', 'Shikimates_and_Phenylpropanoids',
               'Alkaloids')
deduplicated_genus_tree = ape::read.tree(file.path('inputs','prepared_final_smbtree.tre'))
genus_level_data = read.csv(file.path('..','collect_compound_data','outputs','genus_level_pathway_data.csv'))
labelled_tree = get_subset_of_tree_from_genera_in_data(genus_level_data,deduplicated_genus_tree)

for (pathway in np_pathways){
  
  calculate_signal(deduplicated_genus_tree,genus_level_data,paste('mean_identified_as',pathway,sep='_'))
  calculate_signal(deduplicated_genus_tree,genus_level_data,paste('norm_mean_identified_as',pathway,sep='_'))
}



pathway_plot <- function(vars_to_analyse){
  data_with_tree_labels=get_matching_genus_labels(labelled_tree,genus_level_data)[c('label',vars_to_analyse)]
  
  # Convert the dataframe to a data matrix
  data_matrix <- data_with_tree_labels %>%
    select(-label) %>% # Replace "Species_column" with the column name containing species names
    as.matrix()
  rownames(data_matrix) <- data_with_tree_labels$label
  # Associate the data matrix with tree
  circ <- ggtree(labelled_tree, layout = "circular")+geom_tiplab2(size=2, show.legend=FALSE)
  p <- gheatmap(circ, data_matrix, offset=20, width=.2,font.size=1,
                colnames_angle=0, colnames_offset_y = .25) +
    scale_fill_viridis_c(option="A", name="Mean\nIdentified As Pathway")

  output_svg = file.path('outputs', 'pathway_plot.jpg')
  ggplot2::ggsave(output_svg,width=10, height=10,
                  dpi = 300, limitsize=FALSE)
}

pathway_plot(c('mean_identified_as_Polyketides','mean_identified_as_Alkaloids',  'mean_identified_as_Terpenoids','mean_identified_as_Shikimates_and_Phenylpropanoids'))


