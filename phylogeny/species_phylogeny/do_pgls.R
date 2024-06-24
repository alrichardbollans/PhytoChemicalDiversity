library(here)
source(here('helper_functions.R'))

# working_sp_tree = ape::read.tree(file.path('outputs', 'working_species_tree.tre'))
genus_abundance_diversity_data = read.csv(file.path('..','..','diversity_metrics','outputs','genus_level_pathway_diversity_information.csv'), row.names = 1)
genus_distance_diversity_data = read.csv(file.path('..','..','diversity_metrics','outputs','genus_level_distance_diversity_information.csv'), row.names = 1)
all_diversity_data = merge(genus_abundance_diversity_data,genus_distance_diversity_data, by= 'Genus')
# species_data = read.csv(file.path('..','..','collect_compound_data','outputs', 'all_taxa_compound_data.csv'))

phylogenetic_measure_data = read.csv(file.path('outputs', 'phylogenetic_diversities.csv'), row.names = 1)


working_data <- merge(phylogenetic_measure_data,all_diversity_data,by="Genus")
rownames(working_data) <- working_data$Genus

deduplicated_genus_tree = ape::read.tree(genus_tree_path)
working_genus_tree  = subset_tree(deduplicated_genus_tree, working_data$Genus)


plot(working_genus_tree)
geiger::name.check(working_genus_tree, working_data)




working_genus_tree$node.label <- NULL

comp.data<-caper::comparative.data(working_genus_tree, working_data, names.col="Genus", warn.dropped=TRUE)
for(metric in metrics){
  formul = as.formula(sprintf("%s ~ %s", metric, "number_of_species_in_data_and_tree"))

  pglsmodel<-caper::pgls(formul, data=comp.data)
  summary(pglsmodel)
  coef(pglsmodel)
  plot(working_data[, c("number_of_species_in_data_and_tree", metric)])
  abline(a = coef(pglsmodel)[1], b = coef(pglsmodel)[2])
}

for(metric in metrics){
  formul = as.formula(sprintf("%s ~ %s", metric, "phylogenetic_diversity"))
  
  pglsmodel<-caper::pgls(formul, data=comp.data)
  print(summary(pglsmodel))
  coef(pglsmodel)
  plot(working_data[, c("phylogenetic_diversity", metric)])
  abline(a = coef(pglsmodel)[1], b = coef(pglsmodel)[2])
}

