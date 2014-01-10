#heatmaps for illustration
library(pheatmap)
library(entropy)
m = matrix(c(0,1,0,1,1,0,0,1,0), ncol = 3)

pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE)


tmp = features[1:4,c(5:8)]
pheatmap(as.matrix(tmp), cluster_rows = FALSE, cluster_cols = FALSE)


#small graph heatmap
small_graph = read.table("C:\\Users\\v-anleon\\Desktop\\Tartu_University\\Algorithmics2013\\project\\small_graph.txt", sep = ';')
small_graph = as.matrix(small_graph)
pheatmap()

adjacency = matrix(c(0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,
                     0,0,0,0,1,0,0,0,0,0,0,1,0,0), nrow = 7)
pheatmap(adjacency,cluster_rows = T, cluster_cols = TRUE)

#entropy
freqs = apply(adjacency,1,table)
full_graph_entropy = sum(apply(freqs,2,entropy, unit = "log2"))
average_entropy = full_graph_entropy/7

##node1 out
adjacency_updt = adjacency[-c(1),]
adjacency_updt[1,1]=2
#adjacency_updt[2,2]=2
freqs = apply(adjacency_updt,1,table)

entr = sapply(freqs,entropy, unit = "log2")
full_graph_entropy = sum(entr)
average_entropy = full_graph_entropy/6



adjacency_updt2 = adjacency_updt[-c(4),]
adjacency_updt2[3,5]=3
adjacency_updt2[4,5]=3
adjacency_updt2[5,5]=3
freqs = apply(adjacency_updt2,1,table)
full_graph_entropy = sum(apply(freqs,2,entropy, unit = "log2"))
average_entropy = full_graph_entropy/5

#regroup in GraphScope file
##current regrouping:
### segment 1: 4,6,7
### segment 2: 2,1,3
### segment 3: 5
partitioning = adjacency
freqs = apply(partitioning,1,table)
partitioning_entropy = sum(apply(freqs,2,entropy, unit = "log2"))


  
  
