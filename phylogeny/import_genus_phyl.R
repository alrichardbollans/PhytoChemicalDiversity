library(here)



### Importing

cleaned_paftol_tree = ggtree::read.tree(file.path('inputs','genus_apoc_log_rub_tree.tree'))

### Some smb methods
smb_tree = ggtree::read.tree(file.path('inputs','v0.1','ALLMB.tre'))
get_smb_genus_name_from_tree <-function(given_tree_name){
  new = stringr::str_split(given_tree_name,'_')[[1]][1]
  return(new)
}

clean_smb_tree_names <- function(tree){
  tree$tip.label = as.character(lapply(tree$tip.label,get_smb_genus_name_from_tree))
  return(tree)
}

remove_duplicated_tips <- function(tree){
  # Check for repeated tips
  duplicated_tips = which(duplicated(tree$tip.label))
  if (length(duplicated_tips)>0){
    # Check for repeated tips
    print('duplicated tips will be removed, preserving a single instance')
    
    tree = ape::drop.tip(tree,duplicated_tips)
  }
  
  return(tree)
}

cleaned_smb_tree = clean_smb_tree_names(smb_tree)
cleaned_smb_tree = remove_duplicated_tips(cleaned_smb_tree)
print(length(cleaned_smb_tree$tip.label))
