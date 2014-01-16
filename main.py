from igraph import *
from Queue import PriorityQueue
from write_graph import *
from graphScope import *
import time


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
readEdges("edgelists/edgelist_noise_1pct.txt")
graphSegments = []

while(not edges.empty()):
    edgesForTimestamp = []
    timestamp = edges.queue[0].timestamp
    while(not edges.empty() and edges.queue[0].timestamp == timestamp):
        edge = edges.get()
        edgesForTimestamp.append((edge.v1, edge.v2))
     
    g = Graph.Bipartite(types.values(), edgesForTimestamp, directed=False)
    g.vs["label"] = range(len(types))
    graphSegments.append(g)



timestamp = 1
segments = []
for g in graphSegments:
    print "timestamp", timestamp, "\n"
    result = partitionGraph(g, 5)
    nodesInPartS = result[0]
    nodesInPartD = result[1]
    part = result[2]
    print "\n"
    timestamp += 1
    segments.append(GraphSegment(g, result[0], result[1], result[2]))
        # # write files

    
j = 0
while j < len(segments) - 1:
    resultSegment = graphScope(segments[j], segments[j+1])
    segments[j:j+2] = resultSegment
    if len(resultSegment) == 2:
        j+=1
     



  
