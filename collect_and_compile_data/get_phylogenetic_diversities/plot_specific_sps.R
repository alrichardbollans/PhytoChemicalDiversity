library(here)
library(ggtree)
library(ggplot2) # for the scale_color_manual function

library(dplyr) # for filtering

source(here('import_and_check_tree.R'))
labelled_tree = tree_with_tips_in_species_data

tips_to_highlight1=c('Scleromitrion_verticillatum', 'Morinda_citrifolia')
tips_to_highlight2 = c('Vincetoxicum_indicum', 'Leptadenia_reticulata', 'Rauvolfia_serpentina', 'Calotropis_gigantea', 'Leptopetalum_biflorum', 'Guettarda_speciosa', 'Alstonia_scholaris', 'Ochrosia_oppositifolia', 'Morinda_citrifolia', 'Pavetta_indica', 'Oldenlandia_corymbosa', 'Oldenlandia_umbellata')
tips_to_highlight3 = c('Cynanchum_viminale', 'Carissa_spinarum')
tips_to_highlight4 = c('Cynanchum_viminale', 'Guettarda_speciosa', 'Carissa_spinarum') #ALD
tips_to_highlight5 = c('Centaurium_erythraea', 'Galium_palustre', 'Galium_aparine', 'Sherardia_arvensis', 'Vinca_difformis')#AZO
#COR
tips_to_highlight6 = c('Cruciata_pedemontana', 'Centaurium_erythraea', 'Cruciata_glabra', 'Galium_odoratum', 'Cruciata_laevipes', 'Schenkia_spicata', 'Blackstonia_perfoliata', 'Vincetoxicum_hirundinaria', 'Centaurium_pulchellum', 'Rubia_peregrina', 'Galium_mollugo', 'Galium_palustre', 'Galium_album', 'Galium_spurium', 'Galium_tricornutum', 'Galium_aparine', 'Nerium_oleander', 'Sherardia_arvensis', 'Crucianella_maritima', 'Vinca_minor', 'Vinca_difformis', 'Gentiana_lutea', 'Gentiana_asclepiadea')
tips_to_highlight7 = c('Cynanchum_viminale', 'Guettarda_speciosa', 'Cerbera_manghas', 'Tabernaemontana_coffeoides', 'Ochrosia_oppositifolia', 'Carissa_spinarum') #SEY

get_highlight_nodes_for_tips <- function(tips_to_highlight, mrca_node, tree_data){
  # Find the node paths from each tip to the root
  # Get node IDs for the tips you want to highlight
  highlight_nodes <- tree_data %>%
    filter(label %in% tips_to_highlight)%>%
    filter(isTip) %>%
    pull(node)
  
  # Create a data frame to store the nodes along the paths to the root
  
  diff = c(1)
  while (length(diff)>0) {
    node_data = tree_data %>% 
      filter(node %in% highlight_nodes)
    new_nodes = node_data$parent
    new_nodes = new_nodes[new_nodes != mrca_node]
    diff = setdiff(new_nodes, highlight_nodes)
    highlight_nodes = c(new_nodes, highlight_nodes)
  }
  return(highlight_nodes)
}

species_family_data = read.csv(file.path('..','collect_compound_data','outputs', 'species_in_study.csv'))[c('accepted_species', 'accepted_family')]
colnames(species_family_data) <- c('Species','Family')
family_data = get_matching_species_labels(labelled_tree,species_family_data)[c('label','Family')]
# Convert the dataframe to a data matrix
family_data_matrix <- family_data %>%
  select(-label) %>% # Replace "Species_column" with the column name containing species names
  as.matrix()
rownames(family_data_matrix) <- family_data$label

get_plots_for_tips <- function(tips_to_highlight){
  line_size=0.4
  mrca_node <- ape::getMRCA(labelled_tree, tips_to_highlight)

  
  circ <- ggtree(labelled_tree, layout = "circular", size = line_size)
  
  
  # Convert the tree to a data frame that we can filter
  tr_data <- circ$data
  highlight_nodes = get_highlight_nodes_for_tips(tips_to_highlight,mrca_node,tr_data)
  circ <- circ + 
    geom_tree(aes(color = as.factor(node %in% highlight_nodes)), size = 3*line_size, show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE"='#00000000'))+ geom_point2(aes(subset=node==mrca_node), color='red', size=3)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))#+            # Set colors
    # geom_tiplab(data = tr_data %>% filter(label %in% tips_to_highlight),         # Add labels for specific tips
    #             aes(label = label),
    #             size = 3,
    #             hjust = -0.4)
  
  family_offset=0
  family_width =.1
  
  
  # p1 <- gheatmap(circ, family_data_matrix, offset=family_offset, width=family_width,colnames=TRUE, font.size=0) +
  #   scale_fill_viridis_d(option="D", name="Family")+
  #   theme(legend.position = 'none')
  
  return(circ)
}

circ1 = get_plots_for_tips(tips_to_highlight1)
circ2 = get_plots_for_tips(tips_to_highlight2)
circ3 = get_plots_for_tips(tips_to_highlight3)
circ4 = get_plots_for_tips(tips_to_highlight4)
circ5 = get_plots_for_tips(tips_to_highlight5)
circ6 = get_plots_for_tips(tips_to_highlight6)
circ7 = get_plots_for_tips(tips_to_highlight7)
# plot(circ1)

out = cowplot::plot_grid(circ1,NULL, circ2,NULL,circ3,circ4,NULL,circ5,NULL,circ6,circ7,nrow=3,rel_widths = c(1, 0, 1,0,1), labels = c('KZN', 'LDV', 'ROD', 'ALD','AZO','COR','SEY'),label_size=24, align = "hv")
cowplot::save_plot(file.path('outputs', 'highlight_species.jpg'), out,base_asp=1.2,base_height=8)
