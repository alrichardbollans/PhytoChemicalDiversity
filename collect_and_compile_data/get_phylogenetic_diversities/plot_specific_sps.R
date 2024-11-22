library(here)
library(ggtree)
library(ggplot2) # for the scale_color_manual function

library(dplyr) # for filtering

source(here('import_and_check_tree.R'))

tips_to_highlight=c('Scleromitrion_verticillatum', 'Morinda_citrifolia')
labelled_tree = tree_with_tips_in_species_data
# BiocManager::install("ggtreeExtra")
# library(ggtreeExtra) # Needed for fan layout
# Following https://yulab-smu.top/treedata-book/chapter7.html
mrca_node <- ape::getMRCA(labelled_tree, tips_to_highlight)
# Plot the tree
# Create a circular tree plot without labels
circ <- ggtree(labelled_tree, layout = "circular", size = 0.2)+ geom_point2(aes(subset=node==mrca_node), color='red', size=3)
# grp <- list(KZN     = tips_to_highlight)
# groupOTU(circ, grp, 'KZN') + aes(color=KZN) +
#   theme(legend.position="right")


# Convert the tree to a data frame that we can filter
tree_data <- circ$data

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

# Plot the tree with highlighted branches
circ <- circ + 
  geom_tree(aes(color = as.factor(node %in% highlight_nodes)), size = 0.2, show.legend = FALSE) + # Color branches
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))# +            # Set colors
  # geom_tiplab(data = tree_data %>% filter(label %in% tips_to_highlight),         # Add labels for specific tips
  #             aes(label = label), 
  #             size = 3,
  #             hjust = -0.4)

family_offset=0
family_width =.1

species_family_data = read.csv(file.path('..','collect_compound_data','outputs', 'species_in_study.csv'))[c('accepted_species', 'accepted_family')]
colnames(species_family_data) <- c('Species','Family')
family_data = get_matching_species_labels(labelled_tree,species_family_data)[c('label','Family')]
# Convert the dataframe to a data matrix
family_data_matrix <- family_data %>%
  select(-label) %>% # Replace "Species_column" with the column name containing species names
  as.matrix()
rownames(family_data_matrix) <- family_data$label

p1 <- gheatmap(circ, family_data_matrix, offset=family_offset, width=family_width,colnames=TRUE, font.size=0) +
  scale_fill_viridis_d(option="D", name="Family")


output_svg = file.path('outputs', 'highlight_species.jpg')
plotwidth = 9
plotheight = 7
ggplot2::ggsave(output_svg,width=plotwidth, height=plotheight,
                dpi = 600, limitsize=FALSE)
