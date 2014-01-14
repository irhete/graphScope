from igraph import *
from Queue import PriorityQueue
from sets import Set
from math import log
import numpy as np


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
            entropy += (float(arrP[i]) / totalInstancesP) * log(float(arrQ[i]) / totalInstancesQ, 2)
    return -entropy

def totalCost(nodesInPartS, nodesInPartD, part, sourceNodes):
    # source and destination are arrays of frequencies for partitions
    partitionToPartition = np.zeros(shape=(len(nodesInPartS), len(nodesInPartD)))
    for nodeIndex in range(len(sourceNodes)):
        for neighbor in g.neighbors(sourceNodes[nodeIndex]):
            partitionToPartition[part[sourceNodes[nodeIndex]]][part[neighbor]] += 1
    possibleEdgesFromPartitions = np.array([len(elem) for elem in nodesInPartS]) * destNodesCount
    existingEdgesFromPartitions = np.apply_along_axis( sum, axis=1, arr=partitionToPartition)
    complementaryColumnForPartitions = possibleEdgesFromPartitions - existingEdgesFromPartitions
    partitionToPartition = np.append(partitionToPartition, np.transpose(np.matrix(complementaryColumnForPartitions)), axis = 1)

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
        part[sourceNodes[sourceIndex]] = bestPartition
        if len(nodesInPartS[previousPartition]) == 0:
            del nodesInPartS[previousPartition]
            for sourceNode in sourceNodes:
                if part[sourceNode] > previousPartition:
                    part[sourceNode] -= 1;
    

def searchKL(nodesInPartS, nodesInPartD, part, sourceNodes):
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
                print "after split:", nodesInPartS
                reGroup(nodesInPartS, nodesInPartD, part, sourceNodes)
                print "after update:", nodesInPartS
                
        merged = True
        while(merged):
            merged = False
            currentTotalCost = totalCost(nodesInPartS, nodesInPartD, part, sourceNodes)
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
                    if newTotalCost < currentTotalCost:
                        merged = True
                        nodesInPartS = newNodesInPartS
                        part = newPart
                        currentTotalCost = newTotalCost
                        k -= 1
                    else:
                        j += 1
                i += 1
        print "after merge:", nodesInPartS
    return [nodesInPartS, part]


# start!
types = {}
edges = PriorityQueue()
readEdges("C:/Users/Irene/Desktop/KOOL/Text Algorithms/Project/small_graph.txt")
edgesForTimestamp = []
timestamp = edges.queue[0].timestamp
while(not edges.empty() and edges.queue[0].timestamp == timestamp):
    edge = edges.get()
    edgesForTimestamp.append((edge.v1, edge.v2))

g = Graph.Bipartite(types.values(), edgesForTimestamp, directed=True)
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
result = searchKL(nodesInPartS, nodesInPartD, part, sourceNodes)
nodesInPartS = result[0]
part = result[1]
nodesInPartD = [destNodes[:]]
for destNode in destNodes:
    part[destNode] = 0
print "final: \n source:", nodesInPartS, "\n dest:", nodesInPartD, "\n partitioning", part, "\n"

# searchKL for dest nodes    
result = searchKL(nodesInPartD, nodesInPartS, part, destNodes)
nodesInPartD = result[0]
part = result[1]
print "final: \n source:", nodesInPartS, "\n dest:", nodesInPartD, "\n partitioning", part
# 
# result = searchKL(nodesInPartS, nodesInPartD, part, sourceNodes)
# nodesInPartS = result[0]
# part = result[1]
# print "final:", nodesInPartS, nodesInPartD, part



##
##layout = g.layout("kk")
##plot(g, layout = layout)
  
