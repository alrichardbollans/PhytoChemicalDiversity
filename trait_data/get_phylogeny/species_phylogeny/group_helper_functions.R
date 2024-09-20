get_species_in_tree_from_group <- function(tree, group, species_group_info_){
    x = species_group_info_[species_group_info_$Assigned_group == group, ]
    sps = x$accepted_species

    species_list = c()
    for(sp in sps){
    new_name = gsub(" ", "_",sp)
    species_list = c(species_list, new_name)
    }

  return(species_list)
}

check_singleton_group <- function(tree, group, species_group_info_){
  # Checks is a group has only one species in a tree
  group_species_in_tree = get_species_in_tree_from_group(tree,group, species_group_info_)
  if(length(group_species_in_tree)==1){
    # cat("Not more than one species for", genus, "\n")
    return(TRUE)
  }else{
    return(FALSE)
  }
}

check_group_in_tree <- function(tree, group, species_group_info_){
  # Checks is a genus has only one species in a tree
  group_species_in_tree = get_species_in_tree_from_group(tree,group, species_group_info_)
  if(length(group_species_in_tree)>0){
    # cat("Not more than one species for", genus, "\n")
    return(TRUE)
  }else{
    return(FALSE)
  }
}