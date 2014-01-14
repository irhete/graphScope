#libraries
library(igraph)
library(FNN)
library(seriation)

#reading data in format of edge list
initial_data = read.table("C:\\Users\\v-anleon\\Desktop\\Tartu_University\\Algorithmics2013\\project\\tmp.txt", sep =';')
names(initial_data) = c("vertices1","vertices2","timestamp")
graph = graph.data.frame(initial_data, directed=FALSE)
adj = as.matrix(get.adjacency(graph))
plot(graph)

#density calculation
density = graph.density(graph, loops = FALSE)
sbgrph = induced.subgraph(graph, c(1,2,3))
density = graph.density(sbgrph, loops = FALSE)
density

plot(graph)
plot(sbgrph)

#define number of vertices in a graph
vertices = vcount(graph)
combinations = combn(vertices, 3,simplify = TRUE)
sbgrph = induced.subgraph(graph, combinations[,1])
density = graph.density(sbgrph, loops = FALSE)
density

###
k = seriate(adj)

pimage(adj, main = "Random", colorkey=TRUE)
pimage(adj, k, main = "Reordered", colorkey=TRUE)

get_order(k,1)
get_order(k,2)


entropy(c(1,3))