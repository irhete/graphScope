from igraph import *
from Queue import PriorityQueue
from write_graph import *
from graphScope import *


class Edge:
    def __init__(self, v1, v2, timestamp):
        self.v1 = v1
        self.v2 = v2
        self.timestamp = timestamp

def readEdges(fileName):
    f = open(fileName, 'r')
    for line in f:
        pieces = line.strip().split(';')
        edges.put(Edge(int(pieces[0]),int(pieces[1]),pieces[2]),pieces[2])
        types[int(pieces[0])] = 0
        types[int(pieces[1])] = 1
    f.close()





# start!
types = {}
edges = PriorityQueue()
readEdges("edgelists/tmp.txt")
edgesForTimestamp = []
timestamp = edges.queue[0].timestamp
while(not edges.empty() and edges.queue[0].timestamp == timestamp):
    edge = edges.get()
    edgesForTimestamp.append((edge.v1, edge.v2))
 
g = Graph.Bipartite(types.values(), edgesForTimestamp, directed=False)
g.vs["label"] = range(len(types))



totalEdges = len(g.es)
totalNodes = len(g.vs)
destNodesCount = sum(g.vs["type"]*1)
sourceNodesCount = totalNodes - destNodesCount
 
 
part = [0] * totalNodes
nodesInPartS = []
nodesInPartD = []
sourceNodes = [i for i, x in enumerate(g.vs["type"]) if x == False]
destNodes = [i for i, x in enumerate(g.vs["type"]) if x == True]
nodesInPartS.append(sourceNodes[:])
nodesInPartD.append(destNodes[:])
part = [0] * totalNodes
# nodesInPartS = [[1,3],[0],[2]]
 
#initialize destPartition so that each node is in its own partition
destNodePartition = 0
nodesInPartD = []
for destNode in destNodes:
    nodesInPartD.append([destNode])
    part[destNode] = destNodePartition
    destNodePartition += 1
 
 
# searchKL for source nodes

for i in range(3):
    result = searchKL(g, nodesInPartS, nodesInPartD, part, sourceNodes)
    nodesInPartS = result[0]
    part = result[1]
    if (i == 0):
        nodesInPartD = [destNodes[:]]
        for destNode in destNodes:
            part[destNode] = 0
    print "iteration:", i
    print "source final: \n source:", nodesInPartS, "\n dest:", nodesInPartD, "\n partitioning", part
  
    # searchKL for dest nodes    
    result = searchKL(g, nodesInPartD, nodesInPartS, part, destNodes)
    nodesInPartD = result[0]
    part = result[1]
    print "dest final: \n source:", nodesInPartS, "\n dest:", nodesInPartD, "\n partitioning", part


adj = g.get_adjacency()
# writeMatrixToPnms(adj, 'initial_matrix.pnm', 'partitioned_matrix.pnm', sourceNodes, destNodes, nodesInPartS, nodesInPartD);
writeInitialGraphToFile(adj, sourceNodes, destNodes, '1.txt')
writePartitionedGraphToFile(nodesInPartS, nodesInPartD, part, adj, sourceNodes, destNodes, '2.txt')


# layout = g.layout("kk")
# plot(g, layout = layout)
  
