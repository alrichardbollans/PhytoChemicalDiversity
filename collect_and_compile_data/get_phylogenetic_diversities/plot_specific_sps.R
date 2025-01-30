library(here)
library(ggtree)
library(ggplot2) # for the scale_color_manual function

library(dplyr) # for filtering

source(here('import_and_check_tree.R'))
labelled_tree = tree_with_tips_in_species_data


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

get_plots_for_tips <- function(tips_to_highlight, highlight_colour){
  line_size=0.2
  label_size = 24
  mrca_node <- ape::getMRCA(labelled_tree, tips_to_highlight)
  
  
  circ <- ggtree(labelled_tree, layout = "circular", size = line_size)
  
  
  # Convert the tree to a data frame that we can filter
  tr_data <- circ$data
  highlight_nodes = get_highlight_nodes_for_tips(tips_to_highlight,mrca_node,tr_data)
  circ <- circ + 
    geom_tree(aes(color = as.factor(node %in% highlight_nodes)), size = 6*line_size, show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = highlight_colour, "FALSE"='#00000000'))+ geom_point2(aes(subset=node==mrca_node), color=highlight_colour, size=3)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))#+
    # geom_text(aes(0,0,label='N/A'),size=label_size,x = -Inf, y = Inf, hjust = 0, vjust = 1)
  # Set colors
    # geom_tiplab(data = tr_data %>% filter(label %in% tips_to_highlight),         # Add labels for specific tips
    #             aes(label = label),
    #             size = 3,
    #             hjust = -0.4)
  
  # family_offset=0
  # family_width =.1
  
  
  # p1 <- gheatmap(circ, family_data_matrix, offset=family_offset, width=family_width,colnames=TRUE, font.size=0) +
  #   scale_fill_viridis_d(option="D", name="Family")+
  #   theme(legend.position = 'none')
  
  return(circ)
}

get_species_by_group <- function(group_name) {
  # Read the CSV file
  data <- read.csv(file.path('..','get_diversity_metrics','outputs','group_info',"native_regions.csv"))
  
  # Filter rows where Assigned_group matches the given string
  filtered_data <- data[data$Assigned_group == group_name, ]
  
  # Extract the accepted_species values
  species_list <- filtered_data$accepted_species
  
  out_list = as.character(lapply(species_list, add_underscores_to_name))
  
  # Return the list of accepted_species
  return(out_list)
}


# High outliers in for APWD
ALD_tips = get_species_by_group('ALD')
AZO_tips = get_species_by_group('AZO')
COR_tips = get_species_by_group('COR')
KZN_tips = get_species_by_group('KZN')
LDV_tips = get_species_by_group('LDV')
ROD_tips = get_species_by_group('ROD')
SEY_tips = get_species_by_group('SEY')

circ1 = get_plots_for_tips(ALD_tips,"red")
circ2 = get_plots_for_tips(AZO_tips,"red")
circ3 = get_plots_for_tips(COR_tips,"red")
circ4 = get_plots_for_tips(KZN_tips,"red")
circ5 = get_plots_for_tips(LDV_tips,"red")
circ6 = get_plots_for_tips(ROD_tips,"red")
circ7 = get_plots_for_tips(SEY_tips,"red")
plot(circ1)
# Low outliers in for APWD
#CLS
CLS_tips = get_species_by_group('CLS')
#NZN
NZN_tips = get_species_by_group('NZN')
#AGS
AGS_tips = get_species_by_group('AGS')
#MGO
MSO_tips = get_species_by_group('MSO')

circ8 = get_plots_for_tips(AGS_tips,"blue")
circ9 = get_plots_for_tips(CLS_tips,"blue")
circ10 = get_plots_for_tips(MSO_tips,"blue")
circ11 = get_plots_for_tips(NZN_tips,"blue")
plot(circ10)

# Create a spacer plot
spacer <- ggplot() + theme_void()
out = cowplot::plot_grid(circ1, circ2,circ3,circ4,circ5,circ6,circ7,circ8,circ9,circ10,circ11,nrow=4, labels = c('ALD','AZO','COR','KZN', 'LDV', 'ROD','SEY', 'AGS', 'CLS','MSO', 'NZN'),label_size=24)
cowplot::save_plot(file.path('outputs', 'highlight_species.jpg'), out,base_asp=1.2,base_height=10)
