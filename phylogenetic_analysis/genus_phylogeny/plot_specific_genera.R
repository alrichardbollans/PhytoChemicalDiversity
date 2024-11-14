

library(here)
library(ggtree)
library(ggplot2) # for the scale_color_manual function

library(dplyr) # for filtering
source(here('helper_functions.R'))

tips_to_highlight=c('Scleromitrion', 'Morinda')
deduplicated_genus_tree = ape::read.tree(genus_tree_path)
genus_level_data = read.csv(file.path('..','..','collect_and_compile_data','get_diversity_metrics', 'outputs', 'group_data', 'Genus_transformed.csv'))
names(genus_level_data)[names(genus_level_data) == 'Assigned_group'] <- 'Genus'
labelled_tree = get_subset_of_tree_from_genera_in_data(genus_level_data,deduplicated_genus_tree)

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
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +            # Set colors
  geom_tiplab(data = tree_data %>% filter(label %in% tips_to_highlight),         # Add labels for specific tips
              aes(label = label), 
              size = 2,
              hjust = 0)
# circ

family_offset=10
family_width =.1

genus_family_data = read.csv(file.path('inputs', 'genus_family_list.csv'))[c('taxon_name', 'family')]
colnames(genus_family_data) <- c('Genus','Family')
family_data = get_matching_genus_labels(labelled_tree,genus_family_data)[c('label','Family')]
# Convert the dataframe to a data matrix
family_data_matrix <- family_data %>%
  select(-label) %>% # Replace "Species_column" with the column name containing species names
  as.matrix()
rownames(family_data_matrix) <- family_data$label

p1 <- gheatmap(circ, family_data_matrix, offset=family_offset, width=family_width,colnames=TRUE, font.size=0) +
  scale_fill_viridis_d(option="D", name="Family")


output_svg = file.path('outputs', 'highlight_genera.jpg')
plotwidth = 9
plotheight = 7
ggplot2::ggsave(output_svg,width=plotwidth, height=plotheight,
                dpi = 600, limitsize=FALSE)
