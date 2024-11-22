library(here)
source(here('helper_functions.R'))
source(here('group_helper_functions.R'))
source(here('import_and_check_tree.R'))

p = ggtree::ggtree(species_gents_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('outputs',
                               'species_gents_tree.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)

p = ggtree::ggtree(tree_with_tips_in_species_data,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('outputs',
                               'tree_with_tips_in_species_data.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)