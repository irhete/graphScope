class GraphSegment:
    adjacencyString = ""
    n = 0
    m = 0
    cost = 0
    def __init__(self, adjacencyMatrix):
        for row in adjacencyMatrix:
            adjacencyString += row
        n = len(adjacencyMatrix)
        m = len(adjacencyMatrix[0])
        cost = encodingCost(graph)
    def merge(self, graph):
        return 0
    

def encodingCost(segment):
    return 0
def union(graph1, graph2):
    return 0
def searchKL(segment):
    return 0

segments = []

#segment and graph are instances of GraphSegment
def graphScope(segment, graph):
    unionCost = encodingCost(segment)
    if unionCost - segment.cost < graph.cost:
        segment.merge(graph)
        searchKL(segment)
    else:
        segments.append(graph)
        searchKL(graph)

    
    


