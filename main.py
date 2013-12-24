from bitstring import BitArray
from Queue import PriorityQueue
import sys

class Edge:
    def __init__(self, v1, v2):
        self.v1 = v1
	self.v2 = v2

class GraphSegment:
    adjacencyString = BitArray()
    def __init__(self, adjacencyMatrix):
        for row in adjacencyMatrix:
            self.adjacencyString += row
        self.n = len(adjacencyMatrix)
        self.cost = encodingCost(graph)
    def merge(self, graph):
        return 0

def readEdges(fileName):
    f = open(fileName, 'r')
    nodeIdx = 0
    for line in f:
	pieces = line.split(';')
	if (pieces[0] not in nodes):
	    nodes[pieces[0]] = nodeIdx
	    nodeIdx += 1
	if (pieces[1] not in nodes):
	    nodes[pieces[1]] = nodeIdx
	    nodeIdx += 1
        edges.put(Edge(nodes[pieces[0]],nodes[pieces[1]]),pieces[3])
    f.close()

def encodingCost(segment):
    return 0
def union(graph1, graph2):
    return 0
def searchKL(segment):
    return 0


#segment and graph are instances of GraphSegment
def graphScope(segment, graph):
    unionCost = encodingCost(segment)
    if unionCost - segment.cost < graph.cost:
        segment.merge(graph)
        searchKL(segment)
    else:
        segments.append(graph)
        searchKL(graph)

    
segments = []
edges = PriorityQueue()
nodes = {}
readEdges(sys.argv[1])
i = 0
while(not edges.empty()):
    edges.get()
    i+=1
print "Edges:", i
print "Nodes:", len(nodes)	



