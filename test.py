class Graph:
    adjacencyString = ""
    n = 0
    m = 0
    def __init__(self, adjacencyMatrix):
        for row in adjacencyMatrix:
            adjacencyString += row
        n = len(adjacencyMatrix)
        m = len(adjacencyMatrix[0])

class Segment:
    def merge(self, graph):
        return 0
    

def encodingCost(graph, segment):
    return 0

segments = []

def graphScope(segment, segmentCost, graph):
    unionCost = encodingCost(segment, graph)
    graphCost = encodingCost(graph)
    if unionCost - segmentCost < graphCost:
        segment.merge(graph)
        searchKL(segment)
    else:
        newSegment = Segment(graph)
        segments.append(newSegment)
        searchKL(newSegment)

    
    


