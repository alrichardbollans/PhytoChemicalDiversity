library(here)
source(here('helper_functions.R'))
source(here('group_helper_functions.R'))

species_gents_tree = ape::read.tree(species_tree_path)
testit::assert("Tree is ultrametic", ape::is.ultrametric(species_gents_tree))

tags = c('native_region','native_regions_medicinal_species')

for (it in tags){
  filename = paste(it,'.csv', sep="")
  species_group_info = read.csv(file.path('..', 'get_diversity_metrics', 'outputs', 'group_info', filename))
  
  
  species_in_data = species_group_info$accepted_species
  species_data_with_underscores = c()
  for(sp in species_in_data){
    new_name = gsub(" ", "_",sp)
    species_data_with_underscores = c(species_data_with_underscores, new_name)
  }
  
  tree_with_tips_in_species_data = subset_tree(species_gents_tree, species_data_with_underscores)
  
  groups = unique(species_group_info$Assigned_group)
  
  singletons = c()
  measures=c()
  not_in_tree = c()
  for(group in groups){
    is_single = check_singleton_group(tree_with_tips_in_species_data,group, species_group_info)
    if(is_single){
      singletons = c(singletons, group)
      group_df1 = data.frame('Group'=c(group), 'phylogenetic_diversity'=c(NaN),'group_age'=c(NaN), 'number_of_species_in_data_and_tree'= c(1), row.names=c(group))
      measures = rbind(measures, group_df1)
    } else{
      if(check_group_in_tree(tree_with_tips_in_species_data,group, species_group_info)){
        species_in_tree = get_species_in_tree_from_group(tree_with_tips_in_species_data,group, species_group_info)
        calc = calculate_phylogenetic_diversity(tree_with_tips_in_species_data, species_in_tree)
        
        num_sp_in_tree = length(species_in_tree)
        group_df = data.frame('Group'=c(group), 'phylogenetic_diversity'=c(calc$phy_diversity),'group_age'=c(calc$group_age), 'number_of_species_in_data_and_tree'= c(num_sp_in_tree), row.names=c(group))
        measures = rbind(measures, group_df)
      }
      else{
        not_in_tree = c(not_in_tree, group)
      }
      
    }
  }
  # measures$group_age_std = scale(measures$group_age)
  # measures$phylogenetic_diversity_std = scale(measures$phylogenetic_diversity)
  # measures$number_of_species_in_data_and_tree_std = scale(measures$number_of_species_in_data_and_tree)
  write.csv(measures, file.path('outputs','group_data', filename))
}


