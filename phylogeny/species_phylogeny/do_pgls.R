library(here)
source(here('helper_functions.R'))

genus_abundance_diversity_data = read.csv(file.path('..','..','diversity_metrics','outputs','genus_level_pathway_diversity_information.csv'), row.names = 1)
genus_distance_diversity_data = read.csv(file.path('..','..','diversity_metrics','outputs','genus_level_distance_diversity_information.csv'), row.names = 1)
all_diversity_data = merge(genus_abundance_diversity_data,genus_distance_diversity_data, by= 'Genus')

phylogenetic_measure_data = read.csv(file.path('outputs', 'phylogenetic_diversities.csv'), row.names = 1)

working_data <- merge(phylogenetic_measure_data,all_diversity_data,by="Genus")
rownames(working_data) <- working_data$Genus

deduplicated_genus_tree = ape::read.tree(genus_tree_path)
working_genus_tree  = subset_tree(deduplicated_genus_tree, working_data$Genus)


plot(working_genus_tree)
geiger::name.check(working_genus_tree, working_data)


working_genus_tree$node.label <- NULL

comp.data<-caper::comparative.data(working_genus_tree, working_data, names.col="Genus", warn.dropped=TRUE)
# metric = "APWD"
r_df_all_metrics = c()
aic_df_all_metrics = c()
slope_df_all_metrics = c()
p_df_all_metrics = c()
for(metric in metrics){
  
  ## analyse richness
  richness_formul = as.formula(sprintf("%s ~ %s", metric, "number_of_species_in_data_and_tree_std"))
  richness_pglsmodel<-caper::pgls(richness_formul, data=comp.data)
  
  richness_aic = richness_pglsmodel$aic # AIC (Akaike Information Criterion) to determine which model best fits your data, lower the better
  richness_summary = summary(richness_pglsmodel)
  richness_summary
  richness_R2 = richness_summary$r.squared #proportion of the variance in metric explained by features
  richness_p_t_out_df = data.frame(coef(richness_summary)) # coefficients
  richness_f_stat = data.frame(richness_summary$fstatistic)
  richness_slope = coef(richness_pglsmodel)[2]
  richness_p = richness_p_t_out_df$Pr...t..[2]
  
  ## analyse phylogenetic diversity
  div_formul = as.formula(sprintf("%s ~ %s", metric, "phylogenetic_diversity_std"))
  diversity_pglsmodel<-caper::pgls(div_formul, data=comp.data)
  
  diversity_aic = diversity_pglsmodel$aic # AIC (Akaike Information Criterion) to determine which model best fits your data, lower the better
  diversity_summary = summary(diversity_pglsmodel)
  diversity_summary
  diversity_R2 = diversity_summary$r.squared #proportion of the variance in metric explained by features
  diversity_p_t_out_df = data.frame(coef(diversity_summary)) # coefficients
  diversity_f_stat = data.frame(diversity_summary$fstatistic)
  diversity_slope = coef(diversity_pglsmodel)[2]
  diversity_p = diversity_p_t_out_df$Pr...t..[2]
  
  ## analyse both
  both_formul = as.formula(sprintf("%s ~ %s + %s", metric, "phylogenetic_diversity_std", "number_of_species_in_data_and_tree_std"))
  both_pglsmodel<-caper::pgls(both_formul, data=comp.data)
  
  both_aic = both_pglsmodel$aic # AIC (Akaike Information Criterion) to determine which model best fits your data, lower the better
  both_summary = summary(both_pglsmodel)
  both_summary
  both_R2 = both_summary$r.squared #proportion of the variance in metric explained by features
  both_p_t_out_df = data.frame(coef(both_summary)) # coefficients
  both_f_stat = data.frame(both_summary$fstatistic)
  both_diversity_slope = coef(both_pglsmodel)[2]
  both_richness_slope = coef(both_pglsmodel)[3]
  
  both_richness_p = both_p_t_out_df$Pr...t..[3]
  both_diversity_p = both_p_t_out_df$Pr...t..[2]
  
  ### AIC df

  aic_df = data.frame(col=c(richness_aic, diversity_aic, both_aic), row.names=c('Num. Species', 'Phylo diversity', 'Both'))
  names(aic_df)[1] <- sprintf("AIC for %s",metric)
  if(is.null(aic_df_all_metrics)){
    aic_df_all_metrics = aic_df
  }else{
    aic_df_all_metrics = cbind(aic_df_all_metrics, aic_df)
  }
  
  
  ### R2 df
  r_df = data.frame(col=c(richness_R2, diversity_R2, both_R2), row.names=c('Num. Species', 'Phylo diversity', 'Both'))
  names(r_df)[1] <- sprintf("R2 for %s",metric)
  if(is.null(r_df_all_metrics)){
    r_df_all_metrics = r_df
  }else{
    r_df_all_metrics = cbind(r_df_all_metrics, r_df)
  }
  
  
  ### Slope df
  slope_df = data.frame(col=c(richness_slope, diversity_slope, both_diversity_slope, both_richness_slope), row.names=c('Num. Species', 'Phylo diversity', 'Phylo diversity (Both model)', 'Num. Species (Both model)'))
  names(slope_df)[1] <- sprintf("Slope for %s",metric)
  if(is.null(slope_df_all_metrics)){
    slope_df_all_metrics = slope_df
  }else{
    slope_df_all_metrics = cbind(slope_df_all_metrics, slope_df)
  }
  
  ### pvalue df
  p_df = data.frame(col=c(richness_p, diversity_p, both_diversity_p, both_richness_p), row.names=c('Num. Species', 'Phylo diversity', 'Phylo diversity (Both model)', 'Num. Species (Both model)'))
  names(p_df)[1] <- sprintf("P for %s",metric)
  if(is.null(p_df_all_metrics)){
    p_df_all_metrics = p_df
  }else{
    p_df_all_metrics = cbind(p_df_all_metrics, p_df)
  }
  
  richness_p_t_out_df$model = c('just_richness')
  diversity_p_t_out_df$model = c('just_diversity')
  both_p_t_out_df$model = c('both')
  all_summary = rbind(richness_p_t_out_df,diversity_p_t_out_df,both_p_t_out_df)
  
  write.csv(all_summary, file.path('outputs', 'pgls', paste(metric,'.csv', sep='')))
}

write.csv(p_df_all_metrics, file.path('outputs', 'pgls', 'p_values.csv'))

write.csv(r_df_all_metrics, file.path('outputs', 'pgls', 'r2_values.csv'))

write.csv(aic_df_all_metrics, file.path('outputs', 'pgls', 'aic_values.csv'))

write.csv(slope_df_all_metrics, file.path('outputs', 'pgls', 'slope_values.csv'))

