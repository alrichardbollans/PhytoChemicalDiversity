library(here)
source(here('helper_functions.R'))

species_gents_tree = ape::read.tree(species_tree_path)
genus_abundance_diversity_data = read.csv(file.path('..','..','get_diversity_metrics','outputs','genus_level_pathway_diversity_information.csv'), row.names = 1)

species_data = read.csv(file.path('..','..','collect_compound_data','outputs', 'all_species_compound_data.csv'))

species_in_compound_data = species_data$accepted_species

species_in_compound_data_with_underscores = c()
for(sp in species_in_compound_data){
  new_name = gsub(" ", "_",sp)
  species_in_compound_data_with_underscores = c(species_in_compound_data_with_underscores, new_name)
}



tree_with_tips_in_species_data = subset_tree(species_gents_tree, species_in_compound_data_with_underscores)

## Very small number of species that are in compound data but not in tree
species_in_compound_data_that_arent_in_tree = setdiff(species_in_compound_data_with_underscores,tree_with_tips_in_species_data$tip.label)


ape::write.tree(tree_with_tips_in_species_data, file.path('outputs', 'working_species_tree.tre'))

# polyphyletic_genera = c()
# singletons = c()
# measures=c()
# not_in_tree = c()
# for(genus in genus_abundance_diversity_data$Genus){
#   is_single = check_singleton_genus(tree_with_tips_in_species_data,genus)
#   if(is_single){
#     singletons = c(singletons, genus)
#   } else{
#     if(check_genus_in_tree(tree_with_tips_in_species_data, genus)){
#       species_in_tree = get_species_in_tree_from_genus(tree_with_tips_in_species_data,genus)
#       calc = calculate_phylogenetic_diversity(tree_with_tips_in_species_data, species_in_tree)
#
#       genus_species_in_tree = get_species_in_tree_from_genus(tree_with_tips_in_species_data,genus)
#       num_sp_in_tree = length(genus_species_in_tree)
#       genus_df = data.frame('Genus'=c(genus), 'phylogenetic_diversity'=c(calc$phy_diversity),'genus_age'=c(calc$group_age), 'number_of_species_in_data_and_tree'= c(num_sp_in_tree), row.names=c(genus))
#       measures = rbind(measures, genus_df)
#       is_poly = check_polyphyly(tree_with_tips_in_species_data,genus)
#       if(is_poly){
#         polyphyletic_genera= c(polyphyletic_genera,genus)
#       }
#
#     }
#     else{
#       not_in_tree = c(not_in_tree, genus)
#     }
#
#   }
# }
# measures$genus_age_std = scale(measures$genus_age)
# measures$phylogenetic_diversity_std = scale(measures$phylogenetic_diversity)
# measures$number_of_species_in_data_and_tree_std = scale(measures$number_of_species_in_data_and_tree)

# write.csv(measures, file.path('outputs', 'phylogenetic_diversities.csv'))

## Sanity checks
for(g in polyphyletic_genera){
  species_in_tree = get_species_in_tree_from_genus(tree_with_tips_in_species_data,g)
  mrca_node <- ape::getMRCA(tree_with_tips_in_species_data, species_in_tree)
  
  # Extract the clade from the tree
  subtree <- ape::extract.clade(tree_with_tips_in_species_data, node = mrca_node)
  plot(subtree)
}

## Calculation example:
genus = 'Cinchona'
species_in_tree = get_species_in_tree_from_genus(tree_with_tips_in_species_data,genus)
subtree = get_induced_tree_from_species(tree_with_tips_in_species_data, species_in_tree)
plot(subtree)
# branching_times = ape::branching.times(subtree) # Distance of each node to the tips
# 
# branching_times_sum <- sum(ape::branching.times(subtree))
# 
# root_dist = adephylo::distRoot(subtree, species_in_tree, method="patristic") # distances to each tip
# root_dist[1]

total_branch_length = sum(subtree$edge.length) # Faiths measure


## For singleton tips? Ignore
# Find the index of the tip
# sps = get_species_in_tree_from_genus(tree_with_tips_in_species_data,'Aganosma')
# tip_label = sps[1]
# tip_index <- which(tree_with_tips_in_species_data$tip.label == tip_label)
# 
# # The edge matrix in 'ape' has two columns: the first for the parent node and the second for the child node.
# # We want to find the edge where the child node is our tip.
# 
# # Find the edge corresponding to the tip
# tip_edge <- which(tree_with_tips_in_species_data$edge[, 2] == tip_index)
# 
# # Get the edge length
# edge_length <- tree_with_tips_in_species_data$edge.length[tip_edge]
# 
# # Get the parent node (the node that connects to the tip)
# parent_node <- tree_with_tips_in_species_data$edge[tip_edge, 1]
# # Extract the clade from the tree
# subtree <- ape::extract.clade(tree_with_tips_in_species_data, node = parent_node)
# plot(subtree)
