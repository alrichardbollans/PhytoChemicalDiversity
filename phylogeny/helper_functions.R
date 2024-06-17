library(dplyr)


library(ggtree)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggnewscale)

repo_path = Sys.getenv('KEWSCRATCHPATH')
genus_tree_path = file.path(repo_path, 'gentianales_trees','WCVP_12','Smith_and_Brown_ALLMB', 'Genus', 'outputs', 'final_SMB_Gentianales_genus_tree.tre')

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

##### Signal
calculate_signal <- function(tree,genus_data,var_to_analyse){
  genus_data = genus_data[c('Genus',var_to_analyse)]
  
  data_with_tree_labels=get_matching_genus_labels(tree,genus_data)
  data_with_tree_labels[[var_to_analyse]] = as.numeric(data_with_tree_labels[[var_to_analyse]])
  out_dir = file.path('outputs',var_to_analyse)
  dir.create(out_dir)
  labelled_tree = get_subset_of_tree_from_genera_in_data(genus_data,tree)
  reordered_data = reorder_data_frame_like_tree(data_with_tree_labels,labelled_tree)
  
  # Following munkemuller_how_2012, use Pagel's lambda
  signal_lambda <- phytools::phylosig(labelled_tree, reordered_data[[var_to_analyse]], method="lambda", test=TRUE)
  
  # put all measures together in table. D has two p-values
  signal_table <-
    tribble(
      ~metric, ~value, ~pvalue,
      "lambda",  signal_lambda$lambda, signal_lambda$P
    )
  
  # save the table
  write.csv(signal_table, file.path(out_dir,"Genus_phylogenetic_signal_results.csv"))
}


generic_variable_plot <- function(tree,genus_data,var_to_analyse){
  
  genus_data = genus_data[c('Genus',var_to_analyse)]
  
  data_with_tree_labels=get_matching_genus_labels(tree,genus_data)
  data_with_tree_labels[[var_to_analyse]] = as.numeric(data_with_tree_labels[[var_to_analyse]])
  out_dir = file.path('outputs',var_to_analyse)
  dir.create(out_dir)
  
  labelled_tree = get_subset_of_tree_from_genera_in_data(genus_data,tree)
  print('Number of matching tips in labelled tree:')
  print(length(labelled_tree$tip.label))
  tree_to_plot <- full_join(labelled_tree, data_with_tree_labels, by = 'label')
  p = ggtree::ggtree(tree_to_plot,aes(colour=!!sym(var_to_analyse)),layout="circular") + # This was giving an error because of the evaluation of var_to_analyse. Corrected following https://stackoverflow.com/a/53168593/8633026
    ggtree::geom_tiplab2(size=label_size, show.legend=FALSE)+
    scale_color_continuous(low='darkgreen', high='red') #+
  theme(legend.position="right")
  labs(color = var_to_analyse)
  output_svg = file.path(out_dir,'labelled_tree.jpg')
  ggplot2::ggsave(output_svg,width=10, height=5,
                  dpi = dpi, limitsize=FALSE)
}



heatmap_plot <- function(labelled_tree,data_with_tree_labels, legend,outfilename, labelsize = 2.2, plotwidth = 10, plotheight=10, tipnames=FALSE, logtrans=FALSE){
  # Following https://yulab-smu.top/treedata-book/chapter7.html
  genus_family_data = read.csv(file.path('inputs', 'genus_family_list.csv'))[c('taxon_name', 'family')]
  colnames(genus_family_data) <- c('Genus','Family')
  family_data = get_matching_genus_labels(labelled_tree,genus_family_data)[c('label','Family')]
  if(tipnames){
    circ <- ggtree(labelled_tree, layout = "rectangular")+geom_tiplab(size=labelsize, show.legend=FALSE)
    family_offset=10
    variable_offset=15
    family_width =.05
    variable_width = .2
  }else{
    circ <- ggtree(labelled_tree, layout = "rectangular")
    family_offset=1
    variable_offset=10
    family_width =0.1
    variable_width = .4
  }
  
  
  # Convert the dataframe to a data matrix
  family_data_matrix <- family_data %>%
    select(-label) %>% # Replace "Species_column" with the column name containing species names
    as.matrix()
  rownames(family_data_matrix) <- family_data$label
  
  p1 <- gheatmap(circ, family_data_matrix, offset=family_offset, width=family_width,colnames=FALSE) +
    scale_fill_viridis_d(option="D", name="Family")
  
  # Convert the dataframe to a data matrix
  data_matrix <- data_with_tree_labels %>%
    select(-label) %>% # Replace "Species_column" with the column name containing species names
    as.matrix()
  rownames(data_matrix) <- data_with_tree_labels$label
  if (logtrans){
    p2 <- p1 + new_scale_fill()
    gheatmap(p2, data_matrix, offset=variable_offset, width=variable_width,font.size=0,
             colnames_angle=90, colnames_offset_y = 10) +
      scale_fill_viridis_c(option="A", name=legend ,
                           trans = scales::log_trans())
  }else{
    p2 <- p1 + new_scale_fill()
    gheatmap(p2, data_matrix, offset=variable_offset, width=variable_width,font.size=0,
             colnames_angle=90, colnames_offset_y = 10) +
      scale_fill_viridis_c(option="A", name=legend)
  }
  
  
  output_svg = file.path('outputs', outfilename)
  ggplot2::ggsave(output_svg,width=plotwidth, height=plotheight,
                  dpi = 600, limitsize=FALSE)
}
