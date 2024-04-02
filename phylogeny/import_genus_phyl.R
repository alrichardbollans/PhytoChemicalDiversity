library(here)
### Importing
cleaned_paftol_tree = ggtree::read.tree(file.path('inputs','genus_apoc_log_rub_tree.tree'))



get_matching_genus_labels <- function(tree,data){
  # Gets data which appears in tree and appends 'label' column
  # First match by accepted names
  accepted_label_matches <-
    tree$tip.label %>%
    tibble::enframe(name=NULL, value="label")%>%
    rowwise()  %>%
    mutate(Genus=label)%>%
    right_join(
      data,
      by=c("Genus"="Genus")
    )
  
  matching_labels = accepted_label_matches$label
  
  data_with_tree_labels_no_nan = tidyr::drop_na(accepted_label_matches,'label')
  
  return(data_with_tree_labels_no_nan)
}
subset_tree <- function(tree, node_list) {
  drop_list <- tree$tip.label[! tree$tip.label %in% node_list]
  
  return(ape::drop.tip(tree, drop_list))
}

get_subset_of_tree_from_genera_in_data <- function(data, tree){
  
  lab_data = data.frame(data)
  
  # print(lab_data)
  labels = get_matching_genus_labels(tree,lab_data)$label
  
  # drop all tips we haven't found matches for
  f_tree <- subset_tree(tree, labels)
  
  return(f_tree)
}

##### Plotting
library(dplyr)
library(ggtree)
library(ggplot2)
label_size = 2
dpi=600
width=14
height=14

genus_abundance_diversity_data = read.csv(file.path('..','library_info_and_data_import','outputs','genus_level_pathway_diversity_information.csv'))
genus_distance_diversity_data = read.csv(file.path('..','library_info_and_data_import','outputs','genus_level_distance_diversity_information.csv'))

reorder_data_frame_like_tree<-function(df,tree){

  # Get the tip labels from the tree
  tip_labels <- tree$tip.label
  
  # Create a new data frame with all tip labels and their order
  new_data <- data.frame(Genus = tip_labels, row.names = NULL)
  
  # Join the new data frame with the original data
  # This will add rows for missing labels and preserve the order
  merged_data <- left_join(new_data, df, by = "Genus")
  
  if(!setequal(merged_data$Genus,tree$tip.label)){
    stop("Mismatch with data and tree labels.")
  }
  
  return(merged_data)
}

plot_and_calculate_signal <- function(genus_data,var_to_analyse){

  genus_data = genus_data[c('Genus',var_to_analyse)]
  
  data_with_tree_labels=get_matching_genus_labels(cleaned_paftol_tree,genus_data)
  data_with_tree_labels[[var_to_analyse]] = as.numeric(data_with_tree_labels[[var_to_analyse]])
  out_dir = file.path('outputs',var_to_analyse)
  dir.create(out_dir)
  output_svg = file.path(out_dir,'tree.jpg')
  
  tree_to_plot <- full_join(cleaned_paftol_tree, data_with_tree_labels, by = 'label')
  p = ggtree::ggtree(tree_to_plot,layout="circular",
                     mapping=aes(colour=!!sym(var_to_analyse))) +# This was giving an error because of the evaluation of var_to_analyse. Corrected following https://stackoverflow.com/a/53168593/8633026
    geom_tiplab2(size=label_size, show.legend=FALSE)+
    scale_color_continuous(low='darkgreen', high='red') +
    theme(legend.position="right")
    labs(color = var_to_analyse)
  ggplot2::ggsave(output_svg,width=width, height=height,
                  dpi = dpi, limitsize=FALSE)
  
  labelled_tree = get_subset_of_tree_from_genera_in_data(genus_data,cleaned_paftol_tree)
  tree_to_plot <- full_join(labelled_tree, data_with_tree_labels, by = 'label')
  p = ggtree::ggtree(tree_to_plot,aes(colour=!!sym(var_to_analyse)),layout="circular") + # This was giving an error because of the evaluation of var_to_analyse. Corrected following https://stackoverflow.com/a/53168593/8633026
    geom_tiplab2(size=label_size, show.legend=FALSE)+
    scale_color_continuous(low='darkgreen', high='red') #+
    theme(legend.position="right")
    labs(color = var_to_analyse)
  output_svg = file.path(out_dir,'labelled_tree.jpg')
  ggplot2::ggsave(output_svg,width=width, height=height,
                  dpi = dpi, limitsize=FALSE)
  output_svg = file.path(out_dir,'labelled_tree.jpg')
  ggplot2::ggsave(output_svg,width=width, height=height,
                  dpi = dpi, limitsize=FALSE)
  ##### Signal
  
  # calculate Blomberg's K and Pagel's lambda
  
  reordered_data = reorder_data_frame_like_tree(data_with_tree_labels,labelled_tree)

  signal_K <- phytools::phylosig(labelled_tree, reordered_data[[var_to_analyse]], method="K", test=TRUE)
  signal_lambda <- phytools::phylosig(labelled_tree, reordered_data[[var_to_analyse]], method="lambda", test=TRUE)

  # put all measures together in table. D has two p-values
  signal_table <-
    tribble(
      ~metric, ~value, ~pvalue,
      "K",  signal_K$K, signal_K$P,
      "lambda",  signal_lambda$lambda, signal_lambda$P
    )

  # save the table
  write.csv(signal_table, file.path(out_dir,"Genus_phylogenetic_signal_results.csv"))
}

plot_and_calculate_signal(genus_abundance_diversity_data,'bc_shannon')
plot_and_calculate_signal(genus_distance_diversity_data,'N')
plot_and_calculate_signal(genus_distance_diversity_data,'APWD')
plot_and_calculate_signal(genus_abundance_diversity_data,'pielou')
plot_and_calculate_signal(genus_abundance_diversity_data,'simpson')
