repo_path = Sys.getenv('KEWSCRATCHPATH')
species_tree_path = file.path(repo_path, 'gentianales_trees','WCVP_12','Uphy',  'outputs', 'Species', 'Uphylomaker_species_tree.tre')
genus_tree_path = file.path(repo_path, 'gentianales_trees','WCVP_12','Uphy',  'outputs', 'Genus', 'Uphylomaker_genus_tree.tre')

metrics = c('H', 'Hbc', 'G', 'J', 'FAD', 'MFAD', 'APWD')

# Get all species in tree for genus

get_species_in_tree_from_genus <- function(tree, genus){
  species_list <- grep(paste("^",genus,"_", sep=""), tree$tip.label, value = TRUE)
  return(species_list)
}


get_induced_tree_from_species <- function(tree, species_list){
  if(length(species_list)<2){
    cat("Not more than one species for", species_list, "\n")
  } else{
    # mrca_node <- ape::getMRCA(tree, species_list)
    # 
    # # Extract the clade from the tree
    # subtree <- ape::extract.clade(tree, node = mrca_node)
    subtree = ape::keep.tip(tree, species_list)
    return(subtree)
  }
  
}

check_polyphyly <- function(tree, genus){
  # print(genus)
  
  genus_species_in_tree = get_species_in_tree_from_genus(tree,genus)
  if(length(genus_species_in_tree)<2){
    is_polyphyletic <- FALSE
  }else{
    mrca_node <- ape::getMRCA(tree, genus_species_in_tree)

    # Extract the clade from the tree
    subtree <- ape::extract.clade(tree, node = mrca_node)
    # Get the tip labels of the subtree
    subtree_species <- subtree$tip.label
    
    # Check if the subtree species match the species list
    is_polyphyletic <- !all(subtree_species %in% genus_species_in_tree)
    
  }
  
  
  
  if(is_polyphyletic){
    polyphyly_species = setdiff(subtree_species, genus_species_in_tree)
    # cat("Is the group polyphyletic?", is_polyphyletic, "\n")
    # print(genus)
    # print(polyphyly_species)
  }
  
  return(is_polyphyletic)
}

check_singleton_genus <- function(tree, genus){
  # Checks is a genus has only one species in a tree
  genus_species_in_tree = get_species_in_tree_from_genus(tree,genus)
  if(length(genus_species_in_tree)==1){
    # cat("Not more than one species for", genus, "\n")
    return(TRUE)
  }else{
    return(FALSE)
  }
}

check_genus_in_tree <- function(tree, genus){
  # Checks is a genus has only one species in a tree
  genus_species_in_tree = get_species_in_tree_from_genus(tree,genus)
  if(length(genus_species_in_tree)>0){
    # cat("Not more than one species for", genus, "\n")
    return(TRUE)
  }else{
    return(FALSE)
  }
}

calculate_phylogenetic_diversity <- function(tree, species_in_tree){
  subtree = get_induced_tree_from_species(tree, species_in_tree)
  # Calculate the phylogenetic diversity
  # phy_diversity <- sum(ape::branching.times(subtree))
  genus_age <-adephylo::distRoot(subtree, species_in_tree, method="patristic")[1]
  # Faiths measure (Faith 1992)
  phy_diversity <- sum(subtree$edge.length)
  #cat("Phylogenetic Diversity of",genus, ":", phy_diversity, "\n")

  return(list("phy_diversity" = phy_diversity, "group_age"=genus_age))
  
}

subset_tree <- function(tree, sp_list) {
  # Drop tips in tree that aren't in sp_list
  drop_list <- tree$tip.label[! tree$tip.label %in% sp_list]
  
  return(ape::drop.tip(tree, drop_list))
}

add_underscores_to_name <- function(sp_name){
  new_name = gsub(" ", "_",sp_name)
  return(new_name)
}

get_matching_species_labels <- function(tree,data){
  # Gets data which appears in tree and appends 'label' column
  # First match by accepted names
  data['Species'] <- lapply(data['Species'], add_underscores_to_name)
  accepted_label_matches <-
    tree$tip.label %>%
    tibble::enframe(name=NULL, value="label")%>%
    rowwise()  %>%
    mutate(Species=label)%>%
    right_join(
      data,
      by=c("Species"="Species")
    )
  
  matching_labels = accepted_label_matches$label
  
  data_with_tree_labels_no_nan = tidyr::drop_na(accepted_label_matches,'label')
  
  return(data_with_tree_labels_no_nan)
}
