library(here)
source(here('helper_functions.R'))

species_gents_tree = ape::read.tree(species_tree_path)
genus_abundance_diversity_data = read.csv(file.path('..','..','diversity_metrics','outputs','genus_level_pathway_diversity_information.csv'), row.names = 1)

species_data = read.csv(file.path('..','..','collect_compound_data','outputs', 'all_taxa_compound_data.csv'))

species_in_compound_data = species_data$accepted_name

species_in_compound_data_with_underscores = c()
for(sp in species_in_compound_data){
  new_name = gsub(" ", "_",sp)
  species_in_compound_data_with_underscores = c(species_in_compound_data_with_underscores, new_name)
}




tree_with_tips_in_species_data = subset_tree(species_gents_tree, species_in_compound_data_with_underscores)

ape::write.tree(tree_with_tips_in_species_data, file.path('outputs', 'working_species_tree.tre'))

polyphyletic_genera = c()
for(genus in genus_abundance_diversity_data$Genus){
  is_poly = check_polyphyly(tree_with_tips_in_species_data,genus)
  if(is_poly){
    polyphyletic_genera= c(polyphyletic_genera,genus)
  }
}


singletons = c()
measures=c()
for(genus in genus_abundance_diversity_data$Genus){
  is_single = check_singleton_genus(tree_with_tips_in_species_data,genus)
  if(is_single){
    singletons = c(singletons, genus)
  } else{
    if(check_genus_in_tree(tree_with_tips_in_species_data, genus)){
      phy_diversity = calculate_phylogenetic_diversity(tree_with_tips_in_species_data, genus)
      
      genus_species_in_tree = get_species_in_tree_from_genus(tree_with_tips_in_species_data,genus)
      num_sp_in_tree = length(genus_species_in_tree)
      genus_df = data.frame('Genus'=c(genus), 'phylogenetic_diversity'=c(phy_diversity), 'number_of_species_in_data_and_tree'= c(num_sp_in_tree), row.names=c(genus))
      measures = rbind(measures, genus_df)
      
      
    }
    
  }
}

write.csv(measures, file.path('outputs', 'phylogenetic_diversities.csv'))