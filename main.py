from igraph import *
from Queue import PriorityQueue
from sets import Set
from math import log
import numpy as np
from string import *
from networkx.generators.bipartite import bipartite_random_graph


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


##def graphScope(segment, encodingCost0, newGraph):
##    encodingCostUnion
##    encodingCostNewGraph
##    if encodingCostUnion - encodingCost0 < encodingCostNewGraph:
##        segment.append(newGraph)
##        searchKL(segment)
##    else:

def entropy(arr):
    entropy = 0
    totalInstances = sum(arr)
    for instance in arr:
        if instance > 0:
            entropy += (float(instance) / totalInstances) * log(float(instance) / totalInstances, 2)
    return -entropy

def crossEntropy(arrP, arrQ):
    entropy = 0
    totalInstancesP = sum(arrP)
    totalInstancesQ = sum(arrQ)
    for i in range(len(arrP)):
        if arrP[i] > 0:
            if (arrQ[i] == 0):
                print arrP
                print arrQ
            entropy += (float(arrP[i]) / totalInstancesP) * log(float(arrQ[i]) / totalInstancesQ, 2)
    return -entropy

def totalCost(nodesInPartS, nodesInPartD, part, sourceNodes):
    # source and destination are arrays of frequencies for partitions
    partitionToPartition = createPartitionToPartition(nodesInPartS, nodesInPartD, part, sourceNodes)

    source = [len(elem) for elem in nodesInPartS]
    destination = [len(elem) for elem in nodesInPartD]
    m = sum(source)
    n = sum(destination)
    partitionEncodingCost = m * entropy(source) + n * entropy(destination)
    graphEncodingCost = totalEdges
    k = len(source)
    l = len(destination)
    for i in range(k):
        for j in range(l):
            gamma = source[i] * destination[j]
            arr = []
            arr.append(partitionToPartition[i].tolist()[0][j])
            arr.append(gamma - arr[0])
            graphEncodingCost += gamma * entropy(arr)
    return partitionEncodingCost + graphEncodingCost

def averageEntropy(nodesInPartS, nodesInPartD, part, partIndex):
    partEntropy = 0
    for node in nodesInPartS[partIndex]:
        arr = [0] * len(nodesInPartD)
        for neighbor in g.neighbors(node):
            arr[part[neighbor]] += 1
        arr.append(destNodesCount - sum(arr))
        partEntropy += entropy(arr)
    return partEntropy / len(nodesInPartS)

def findPartitionToSplit(nodesInPartS, nodesInPartD, part):
    maxEntropy = 0
    maxPartition = 0
    for partIndex in range(len(nodesInPartS)):
        avgEntropy = averageEntropy(nodesInPartS, nodesInPartD, part, partIndex)
#         print "average entropy for partition", partIndex, "=", avgEntropy
        if avgEntropy > maxEntropy:
            maxEntropy = avgEntropy;
            maxPartition = partIndex
    return maxPartition

def reGroup(nodesInPartS, nodesInPartD, part, sourceNodes):
    sourceNodesCount = len(sourceNodes)
    destNodesCount = totalNodes - len(sourceNodes)
    partitionToPartition = np.zeros(shape=(len(nodesInPartS), len(nodesInPartD)))
    nodeToPartition = np.zeros(shape=(sourceNodesCount, len(nodesInPartD)))
    for node in range(sourceNodesCount):
        for neighbor in g.neighbors(sourceNodes[node]):
            nodeToPartition[node][part[neighbor]] += 1
            partitionToPartition[part[sourceNodes[node]]][part[neighbor]] += 1
    complementaryColumnForNodes = destNodesCount - np.apply_along_axis( sum, axis=1, arr=nodeToPartition)
    nodeToPartition = np.append(nodeToPartition, np.transpose(np.matrix(complementaryColumnForNodes)), axis = 1)
    possibleEdgesFromPartitions = np.array([len(elem) for elem in nodesInPartS]) * destNodesCount
   
    existingEdgesFromPartitions = np.apply_along_axis( sum, axis=1, arr=partitionToPartition)
    complementaryColumnForPartitions = possibleEdgesFromPartitions - existingEdgesFromPartitions
    partitionToPartition = np.append(partitionToPartition, np.transpose(np.matrix(complementaryColumnForPartitions)), axis = 1)

    for sourceIndex in range(sourceNodesCount):
        bestPartitionCost = 1000000
        for partition in range(len(nodesInPartS)):
            if (partition == part[sourceNodes[sourceIndex]]):
                currentCrossEntropy = crossEntropy(nodeToPartition[sourceIndex].tolist()[0], partitionToPartition[partition].tolist()[0])
            else:
                currentCrossEntropy = crossEntropy(nodeToPartition[sourceIndex].tolist()[0], (partitionToPartition[partition] + nodeToPartition[sourceIndex]).tolist()[0])
#             print "node:", sourceNodes[sourceIndex], ", partition:", partition, ", cost:", currentCrossEntropy
            if currentCrossEntropy < bestPartitionCost:
                bestPartition = partition
                bestPartitionCost = currentCrossEntropy
        previousPartition = part[sourceNodes[sourceIndex]]
        nodesInPartS[previousPartition].remove(sourceNodes[sourceIndex])
        nodesInPartS[bestPartition].append(sourceNodes[sourceIndex])
        partitionToPartition[previousPartition] = partitionToPartition[previousPartition] - nodeToPartition[sourceIndex];
        partitionToPartition[bestPartition] = partitionToPartition[bestPartition] + nodeToPartition[sourceIndex];
        part[sourceNodes[sourceIndex]] = bestPartition
        if len(nodesInPartS[previousPartition]) == 0:
            del nodesInPartS[previousPartition]
            partitionToPartition = np.delete(partitionToPartition, (previousPartition), axis=0)
            for sourceNode in sourceNodes:
                if part[sourceNode] > previousPartition:
                    part[sourceNode] -= 1;
    

def searchKL(nodesInPartS, nodesInPartD, part, sourceNodes, verbose = False):
    changed = True
    while (changed):
        changed = False

        encodingCostDecreased = True
        while(encodingCostDecreased):
            encodingCostDecreased = False
            partIndex = findPartitionToSplit(nodesInPartS, nodesInPartD, part)
            if len(nodesInPartS[partIndex]) > 1:
                currentAverageEntropy = averageEntropy(nodesInPartS, nodesInPartD, part, partIndex)
#                 print "current average entropy", currentAverageEntropy
                for s in nodesInPartS[partIndex]:
                    newNodesInPartS = nodesInPartS[:]
                    newNodesInPartS[partIndex].remove(s)
                    newNodesInPartS.append([s])
                    newPart = part[:]
                    newPart[s] = len(nodesInPartS)
                    newAverageEntropy = averageEntropy(newNodesInPartS, nodesInPartD, newPart, partIndex)
#                     print "new average entropy", newAverageEntropy
                    if newAverageEntropy < currentAverageEntropy:
                        nodesInPartS = newNodesInPartS
                        currentAverageEntropy = newAverageEntropy
                        part = newPart
                if (verbose):
                    print "after split:", nodesInPartS
                reGroup(nodesInPartS, nodesInPartD, part, sourceNodes)
                if (verbose):
                    print "after update:", nodesInPartS
                
        merged = True
        while(merged):
            merged = False
            currentTotalCost = totalCost(nodesInPartS, nodesInPartD, part, sourceNodes)
#             print "current total cost:", currentTotalCost
            k = len(nodesInPartS)
            i = 0
            while i < k-1:
                j = i+1
                while j < k:
                    newNodesInPartS = []
                    for partition in range(len(nodesInPartS)):
                        newNodesInPartS.append(nodesInPartS[partition][:])
                    newNodesInPartS[i].extend(newNodesInPartS[j])
                    del newNodesInPartS[j]
                    newPart = part[:]
                    for sourceNode in sourceNodes:
                        if newPart[sourceNode] == j:
                            newPart[sourceNode] = i
                        elif newPart[sourceNode] > j:
                            newPart[sourceNode] -= 1
                    newTotalCost = totalCost(newNodesInPartS, nodesInPartD, newPart, sourceNodes)
#                     print "new total cost:", newTotalCost, "source:", nodesInPartS
                    if newTotalCost <= currentTotalCost:
                        merged = True
                        nodesInPartS = newNodesInPartS
                        part = newPart
                        currentTotalCost = newTotalCost
                        k -= 1
                    else:
                        j += 1
                i += 1
        if (verbose):
            print "after merge:", nodesInPartS
    return [nodesInPartS, part]

def createPartitionToPartition(nodesInPartS, nodesInPartD, part, sourceNodes):
    partitionToPartition = np.zeros(shape=(len(nodesInPartS), len(nodesInPartD)))
    for nodeIndex in range(len(sourceNodes)):
        for neighbor in g.neighbors(sourceNodes[nodeIndex]):
            partitionToPartition[part[sourceNodes[nodeIndex]]][part[neighbor]] += 1
    destNodesCount = totalNodes - len(sourceNodes)
    possibleEdgesFromPartitions = np.array([len(elem) for elem in nodesInPartS]) * destNodesCount
    existingEdgesFromPartitions = np.apply_along_axis( sum, axis=1, arr=partitionToPartition)
    complementaryColumnForPartitions = possibleEdgesFromPartitions - existingEdgesFromPartitions
    partitionToPartition = np.append(partitionToPartition, np.transpose(np.matrix(complementaryColumnForPartitions)), axis = 1)
    return partitionToPartition

def findMatrixReordering(nodesInPart):
    ordering = []
    for partition in nodesInPart:
        for node in partition:
            ordering.append(node)
    ordering = np.array(ordering) - min(ordering)
    return ordering

def writeMatrixToFile(adj, fileName, newLines = False):
    f = open(fileName,'a')
    for row in adj:
        f.write("".join(str(x) for x in row))
        if (newLines):
            f.write("\n")
    if (not newLines):
        f.write("\n")
    f.close()
    
def writeArrayToFile(arr, fileName):
    f = open(fileName,'a')
    f.write("".join(str(x) for x in arr))
    f.write("\n")
    f.close()

def writeNumbersToFile(numbers, fileName):
    f = open(fileName, 'w')
    for number in numbers:
        f.write(str(number) + "\n")
    f.close()

def writePartitionedGraphToFile(nodesInPartS, nodesInPartD, part, adj, fileName):
    writeNumbersToFile([sourceNodesCount, destNodesCount, len(nodesInPartS), len(nodesInPartD)], fileName)
    writeArrayToFile(part, fileName)
    for i in range(len(nodesInPartS)):
        nodesInPartS[i].sort()
        for j in range(len(nodesInPartD)):
            nodesInPartD[j].sort()
            adj2 = [[adj[row][col] for col in nodesInPartD[j]] for row in nodesInPartS[i]]
            writeMatrixToFile(adj2, fileName)
            
def writeInitialGraphToFile(adj, fileName):
    writeNumbersToFile([sourceNodesCount, destNodesCount], fileName)
    adj2 = [[adj[row][col] for col in destNodes] for row in sourceNodes]
    writeMatrixToFile(adj2, fileName)

def writeMatrixToPnm(adj, fileName1, fileName2, sourceNodes, destNodes):
    adj2 = [[adj[row][col] for col in destNodes] for row in sourceNodes]
    sourceOrdering = findMatrixReordering(nodesInPartS)
    destOrdering = findMatrixReordering(nodesInPartD)
    adj3 = [[adj2[row][col] for col in destOrdering] for row in sourceOrdering]

    f = open(fileName1,'w')
    f.write("P1\n")
    f.write(str(len(sourceNodes)) + " " + str(len(destNodes)) + "\n")
    f.close()
    writeMatrixToFile(adj2, fileName1, True)
    
    f = open(fileName2,'w')
    f.write("P1\n")
    f.write(str(len(sourceNodes)) + " " + str(len(destNodes)) + "\n")
    f.close()
    writeMatrixToFile(adj3, fileName2, True)
    

# start!
types = {}
edges = PriorityQueue()
readEdges("edgelist (2).txt")
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
    result = searchKL(nodesInPartS, nodesInPartD, part, sourceNodes)
    nodesInPartS = result[0]
    part = result[1]
    if (i == 0):
        nodesInPartD = [destNodes[:]]
        for destNode in destNodes:
            part[destNode] = 0
    print "iteration:", i
    print "source final: \n source:", nodesInPartS, "\n dest:", nodesInPartD, "\n partitioning", part
  
    # searchKL for dest nodes    
    result = searchKL(nodesInPartD, nodesInPartS, part, destNodes)
    nodesInPartD = result[0]
    part = result[1]
    print "dest final: \n source:", nodesInPartS, "\n dest:", nodesInPartD, "\n partitioning", part


adj = g.get_adjacency()
writeMatrixToPnm(adj, 'initial_matrix.pnm', 'partitioned_matrix.pnm', sourceNodes, destNodes);
# writeInitialGraphToFile(adj, 'initialGraph.txt')
# writePartitionedGraphToFile(nodesInPartS, nodesInPartD, part, adj, 'graph.txt')

# print findMatrixReordering(nodesInPartS)
# print findMatrixReordering(nodesInPartD)

    #print matrix
#     adj = g.get_adjacency()
#     writeMatrixToFile(adj, 'graph.txt')



# layout = g.layout("kk")
# plot(g, layout = layout)
  
