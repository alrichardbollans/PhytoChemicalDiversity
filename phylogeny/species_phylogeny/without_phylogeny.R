library(here)
source(here('helper_functions.R'))

genus_abundance_diversity_data = read.csv(file.path('..','..','diversity_metrics','outputs','genus_level_pathway_diversity_information.csv'), row.names = 1)
genus_distance_diversity_data = read.csv(file.path('..','..','diversity_metrics','outputs','genus_level_distance_diversity_information.csv'), row.names = 1)
all_diversity_data = merge(genus_abundance_diversity_data,genus_distance_diversity_data, by= 'Genus')

phylogenetic_measure_data = read.csv(file.path('outputs', 'phylogenetic_diversities.csv'), row.names = 1)

working_data <- merge(phylogenetic_measure_data,all_diversity_data,by="Genus")
rownames(working_data) <- working_data$Genus

# Clearly neither are normal
hist(working_data$phylogenetic_diversity)
hist(working_data$number_of_species_in_data_and_tree)



phyl_diversity_tests = c()
for(metric in metrics){
  div_test = cor.test(working_data$phylogenetic_diversity, working_data[,metric], method = "kendall")
  metric_df = data.frame('Metric'=c(metric), 'method'=c(div_test$method), 'statistic'=c(div_test$statistic), 'pvalue'= c(div_test$p.value), row.names=c(metric))
  phyl_diversity_tests = rbind(metric_df,phyl_diversity_tests)
}
write.csv(phyl_diversity_tests, file.path('outputs', 'correlations_with_phyl_diversity.csv'))

num_sp_tests = c()
for(metric in metrics){
  div_test = cor.test(working_data$number_of_species_in_data_and_tree, working_data[,metric], method = "kendall")
  metric_df = data.frame('Metric'=c(metric), 'method'=c(div_test$method), 'statistic'=c(div_test$statistic), 'pvalue'= c(div_test$p.value), row.names=c(metric))
  num_sp_tests = rbind(metric_df,num_sp_tests)
}
write.csv(num_sp_tests, file.path('outputs', 'correlations_with_num_species.csv'))
