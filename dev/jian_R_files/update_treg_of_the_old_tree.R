library(yaml)
library(data.tree)

new_tree <- yaml.load_file("dev/tree.yaml") %>% 
  as.Node

new_tree$immune_cell$t_cell$t_CD4$t_helper$RemoveChild("t_reg")

new_tree$immune_cell$t_cell$t_CD4$AddChild("t_reg")

Prune(new_tree$epithelial, function(node) node$name == "epithelial")

Sort(new_tree, function(node) node$name)

# save the new tree as yaml file
new_tree %>% 
  as.list %>% 
  write_yaml("dev/jian_R_files/new_tree.yaml")



# load the saved tree
new_tree <- yaml.load_file("dev/jian_R_files/new_tree.yaml") %>% 
  as.Node

