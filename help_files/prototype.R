#libraries
library(igraph)
library(FNN)
library(seriation)
library(reshape)
library(gridExtra)
library(sna)
setwd("C:/Users/v-anleon/Desktop/Tartu_University/Algorithmics2013/project")

reordering = function(file, source_order = source_ordering, destination_order=destination_ordering){
  edgelist = read.table(file, header = FALSE, sep =';')
  names(edgelist) = c("source_nodes","destin_nodes","timestamp")
  source_order = as.character(source_order)
  destination_order = as.character(destination_order)
  bipartite_adj = as.matrix(table(edgelist$source_nodes, edgelist$destin_nodes))
  initial = bipartite_adj[,]
  reordered = bipartite_adj[source_order, destination_order]
  par(mfrow = c(1,2))
  plot.sociomatrix(initial, main = "initial", drawlines = FALSE)
  abline(h = 1.5, v = 1.5, col ="red", lty = 3)  
  plot.sociomatrix(reordered, main = 'reordered', drawlines = FALSE)
  par(mfrow = c(1,1))
}
pairwise_combination = function(v,p){
  u = c()
  for(i in 1:length(v)){
    for(j in 1:length(p)){
      u = rbind(u, c(v[i],p[j]))    
    }
  }
  return(u)
}
split_fn = function(x, groups){
  groups = groups-1
  g <- factor(round(groups * runif(groups * x)))
  bins <- split(x, g)
  return(bins)
}


source_ordering = c(7, 9, 10, 11, 6, 8)
destination_ordering = c(0, 1, 4, 2, 3, 5)
path = "C:/Users/v-anleon/Desktop/Tartu_University/Algorithmics2013/project/tmp.txt"
reordering(path)

edge_freq = rbind(c(0,12),c(6,0))
k = seriate(edge_freq)
get_order(k,1)
get_order(k,2)

edgelist <- cbind.data.frame(source_nodes = sample(0:(m-1), edges, replace = T), destin_nodes = sample((m):(m+n-1), edges, replace=T))
edgelist$timestamp = 1
edgelist = edgelist[-which(duplicated(edgelist)==TRUE),]

##########################
set.seed(round(runif(1, min=2, max = 200),0))
a = 0:999
b = 1000:1999
groups = 4
sources = split_fn(x=a, groups = 4)
destinations = split_fn(x =b,groups = 4)

edgelist = c()
for(i in 1:groups){
  edgelist = rbind(edgelist, pairwise_combination(sources[[i]],destinations[[i]]))
}
edgelist = as.data.frame(edgelist)
head(edgelist)
edgelist = cbind.data.frame(edgelist,1)

#introducing noise
edges_total = nrow(edgelist)
prcnt_noise = 0.25
#delete 0.5% of edges
deletions = round((prcnt_noise/2)*edges_total,0)
#set.seed(6576)
idx_to_delete = sample(nrow(edgelist), deletions)
edgelist = edgelist[-idx_to_delete,]

#add noise
#random sample of existing source nodes
additions = round((prcnt_noise/2)*edges_total,0)
source_add = sample(a, additions, replace = T)
destination_add = sample(b, additions, replace = T)
added_edges = as.data.frame(cbind(source_add,destination_add,1))
names(added_edges) = names(edgelist)
edgelist = rbind.data.frame(edgelist,added_edges)

edgelist = edgelist[-which(duplicated(edgelist)==TRUE),]

write.table(edgelist, "C:/Users/v-anleon/Desktop/Tartu_University/Algorithmics2013/project_bins/data/edgelist_noise_25pct.txt", col.names= FALSE, row.names = FALSE, sep =';')

#reading in the order found by the algorithm
con  <- file("orderings.txt", open = "r")
sources_order <- readLines(con, n = 1, warn = FALSE)
sources_order = strsplit(sources_order, ", ")[[1]]
sources_order = as.numeric(sources_order)

dest_order <- readLines(con, n = 2, warn = FALSE)
dest_order = strsplit(dest_order, ", ")[[1]]
dest_order = as.numeric(dest_order)

source_ordering = sources_order
destination_ordering = dest_order
path = "edgelist.txt"
reordering(path)


#check
table(edgelist$source_nodes)
table(edgelist$destin_nodes)
g = graph.data.frame(edgelist, directed = FALSE)
plot(g, mark.groups=list(0:9),layout=layout.reingold.tilford)
l <-layout.reingold.tilford(g) 
new.layout = cbind(l[,1],c(rep(0,10),rep(1,nrow(l)-10)))
plot(g, mark.groups=list(0:9),layout=new.layout)
g = graph.data.frame(edgelist, directed = FALSE)