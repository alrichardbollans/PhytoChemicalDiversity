library(here)

### Some smb methods

get_smb_genus_name_from_tree <-function(given_tree_name){
  new = stringr::str_split(given_tree_name,'_')[[1]][1]
  return(new)
}

remove_duplicated_tips <- function(tree){
  # Check for repeated tips
  duplicated_tips = which(duplicated(tree$tip.label))
  if (length(duplicated_tips)>0){
    # Check for repeated tips
    print('duplicated tips will be removed, preserving a single instance')
    print(length(duplicated_tips))
    tree = ape::drop.tip(tree,duplicated_tips)
  }
  
  return(tree)
}

### Importing
## Get gentianales tree
smb_tree_ = ggtree::read.tree(file.path('inputs','v0.1','ALLMB.tre'))
gentianales_tree = ape::extract.clade(smb_tree_,c('Gentianales.rn.d8s.tre'))
ape::write.tree(gentianales_tree, file=file.path('inputs','SMB_ALLMB_Gentianales.tre'))
p = ggtree::ggtree(gentianales_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('inputs','SMB_ALLMB_Gentianales.jpg'),width=20, height=16,
               dpi = 300, limitsize=FALSE)

# Then reimport the tree with names standardised by automatchnames
standardised_smb_tree = ggtree::read.tree(file.path('inputs','standardised_smb_tree.tre'))

testit::assert("Check name standardisation preserves tree tips", length(gentianales_tree$tip.label) == length(standardised_smb_tree$tip.label))


## Find nodes and tips of accepted genera
accepted_genus_list <- readLines(file.path('inputs','genus_list.txt'))
common_genera = intersect(accepted_genus_list,standardised_smb_tree$node.label)

#### Add genus nodes as tips

# Function modified from http://blog.phytools.org/2012/11/adding-single-tip-to-tree.html
bind_node_tip<-function(tree,node_label){
  node_index <-match(node_label, tree$node.label)+length(tree$tip.label)
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=node_label,
            edge.length=0,
            Nnode=1)
  class(tip)<-"phylo"
  obj<-ape::bind.tree(tree,tip,where=node_index)
  return(obj)
}

get_descendant_tips <- function(tree, node_label){
  node_index <-match(node_label, tree$node.label)+length(tree$tip.label)
  tips = geiger::tips(tree, node_index)
  tips = tips[! tips %in% c(node_label)] # Don't include node label in descendants
  return(tips)
}

add_node_tip_and_remove_descendant <- function(tree, node_label){
  if(node_label %in% tree$node.label){
    tips_to_drop = get_descendant_tips(tree,node_label)
    new_tree = bind_node_tip(tree,node_label)
    new_tree = ape::drop.tip(new_tree,tips_to_drop)
    return(new_tree)
  }else{
    print('WARNING: Following node no longer in tree')
    print(node_label)
    
    return(tree)
  }
  
}

node_tip_tree = standardised_smb_tree
for (genus in common_genera) {
  node_tip_tree<-add_node_tip_and_remove_descendant(node_tip_tree,genus)
}
new_common_genera = intersect(accepted_genus_list,node_tip_tree$tip.label)
testit::assert("Check extracting nodes preserves important nodes", length(common_genera) < length(new_common_genera))
plot(node_tip_tree)

non_genus_tips_to_remove = node_tip_tree$tip.label[! node_tip_tree$tip.label %in% c(new_common_genera)]
genus_tree = ape::drop.tip(node_tip_tree, non_genus_tips_to_remove)
new_common_genera1 = intersect(accepted_genus_list,genus_tree$tip.label)
testit::assert("Check extracting nodes preserves important nodes", length(new_common_genera) == length(new_common_genera1))
plot(genus_tree)

deduplicated_genus_tree = remove_duplicated_tips(genus_tree)
testit::assert("Check extracting nodes preserves important nodes", length(new_common_genera1) == length(intersect(accepted_genus_list,deduplicated_genus_tree$tip.label)))
p = ggtree::ggtree(deduplicated_genus_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('inputs','SMB_Gentianales_Genus_tree.jpg'),width=20, height=16,
               dpi = 300, limitsize=FALSE)
ape::write.tree(deduplicated_genus_tree, file=file.path('inputs','prepared_final_smbtree.tre'))

