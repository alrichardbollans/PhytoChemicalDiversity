repo_path = Sys.getenv('KEWSCRATCHPATH')
species_tree_path = file.path(repo_path, 'gentianales_trees','WCVP_12','Smith_and_Brown_ALLMB', 'Species', 'outputs', 'final_SMB_Gentianales_species_tree.tre')
genus_tree_path = file.path(repo_path, 'gentianales_trees','WCVP_12','Smith_and_Brown_ALLMB', 'Genus', 'outputs', 'final_SMB_Gentianales_genus_tree.tre')

metrics = c('H', 'Hbc', 'G', 'J', 'FAD', 'MFAD', 'APWD', 'N')

# Get all species in tree for genus

get_species_in_tree_from_genus <- function(tree, genus){
  species_list <- grep(paste("^",genus,"_", sep=""), tree$tip.label, value = TRUE)
  return(species_list)
}


get_induced_tree_from_species <- function(tree, species_list){
  if(length(species_list)<2){
    cat("Not more than one species for", species_list, "\n")
  } else{
    mrca_node <- ape::getMRCA(tree, species_list)
    
    # Extract the clade from the tree
    subtree <- ape::extract.clade(tree, node = mrca_node)
    return(subtree)
  }
  
}

check_polyphyly <- function(tree, genus){
  # print(genus)
  
  genus_species_in_tree = get_species_in_tree_from_genus(tree,genus)
  if(length(genus_species_in_tree)<2){
    is_polyphyletic <- FALSE
  }else{
    subtree = get_induced_tree_from_species(tree, genus_species_in_tree)
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

calculate_phylogenetic_diversity <- function(tree, genus){
  species_in_tree = get_species_in_tree_from_genus(tree,genus)
  subtree = get_induced_tree_from_species(tree, species_in_tree)
  
  
  # Calculate the phylogenetic diversity
  # phy_diversity <- sum(ape::branching.times(subtree))
  # phy_diversity <-adephylo::distRoot(subtree, species_in_tree, method="patristic")[1]
  # Faiths measure (Faith 1992)
  phy_diversity <- sum(subtree$edge.length)
  cat("Phylogenetic Diversity of",genus, ":", phy_diversity, "\n")
  
  return(phy_diversity)
  
}

subset_tree <- function(tree, sp_list) {
  # Drop tips in tree that aren't in sp_list
  drop_list <- tree$tip.label[! tree$tip.label %in% sp_list]
  
  return(ape::drop.tip(tree, drop_list))
}
