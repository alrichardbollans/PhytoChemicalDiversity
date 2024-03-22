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


plot_and_calculate_signal <- function(comp_class){
  genus_compound_data = read.csv(file.path('..','library_info_and_data_import','outputs',paste('genus_level_',comp_class,'_information.csv',sep='')))
  
  data_with_tree_labels=get_matching_genus_labels(cleaned_paftol_tree,genus_compound_data)
  out_dir = file.path('outputs',comp_class)
  dir.create(out_dir)
  output_svg = file.path(out_dir,'tree.svg')
  p = ggtree::ggtree(cleaned_paftol_tree,layout="circular",
                     mapping=aes(colour=mean_identified_as_class)) %<+% data_with_tree_labels +
    geom_tiplab2(size=label_size, show.legend=FALSE)+
    scale_color_gradient2(mid = "blue", high = "red")+ # the lows are nans?
    labs(color = comp_class)
  ggplot2::ggsave(output_svg,width=width, height=height,
                  dpi = dpi, limitsize=FALSE)
  
  labelled_tree = get_subset_of_tree_from_genera_in_data(genus_compound_data,cleaned_paftol_tree)
  
  p = ggtree::ggtree(labelled_tree,layout="circular",
                     mapping=aes(colour=mean_identified_as_class)) %<+% data_with_tree_labels +
    geom_tiplab2(size=label_size, show.legend=FALSE)+
    scale_color_gradient2(mid = "blue", high = "red")+ # the lows are nans?
    labs(color = comp_class)
  output_svg = file.path(out_dir,'labelled_tree.svg')
  ggplot2::ggsave(output_svg,width=width, height=height,
                  dpi = dpi, limitsize=FALSE)
  
  ##### Signal
  
  # calculate Blomberg's K and Pagel's lambda
  
  reordered_data = data_with_tree_labels[match(labelled_tree$tip.label, data_with_tree_labels$Genus),] # Ensure same ordering of data, though this actually happens when data_with_tree_labels is produced
  
  if(!setequal(reordered_data$Genus,labelled_tree$tip.label)){
    stop("Mismatch with data and tree labels. Can probably be fixed by passing named list in phytools methods")
  }
  signal_K <- phytools::phylosig(labelled_tree, reordered_data$mean_identified_as_class, method="K", test=TRUE)
  signal_lambda <- phytools::phylosig(labelled_tree, reordered_data$mean_identified_as_class, method="lambda", test=TRUE)
  
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
plot_and_calculate_signal('alkaloid')
plot_and_calculate_signal('terpenoid')
plot_and_calculate_signal('peptide')
plot_and_calculate_signal('pyrrolizidine')
plot_and_calculate_signal('quinoline')
plot_and_calculate_signal('chemical_diversity')
