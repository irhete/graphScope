from math import log
import numpy as np

# def graphScope(segment, encodingCost0, newGraph):
#     encodingCostUnion
#     encodingCostNewGraph
#     if encodingCostUnion - encodingCost0 < encodingCostNewGraph:
#         segment.append(newGraph)
#         searchKL(segment)
#     else:

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

def totalCost(g, nodesInPartS, nodesInPartD, part, sourceNodes):
    # source and destination are arrays of frequencies for partitions
    partitionToPartition = createPartitionToPartition(g, nodesInPartS, nodesInPartD, part, sourceNodes)

    source = [len(elem) for elem in nodesInPartS]
    destination = [len(elem) for elem in nodesInPartD]
    m = sum(source)
    n = sum(destination)
    partitionEncodingCost = m * entropy(source) + n * entropy(destination)
    graphEncodingCost = len(g.es)
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

def averageEntropy(g, nodesInPartS, nodesInPartD, part, partIndex, sourceNodes):
    partEntropy = 0
    destNodesCount = len(g.vs) - len(sourceNodes)
    for node in nodesInPartS[partIndex]:
        arr = [0] * len(nodesInPartD)
        for neighbor in g.neighbors(node):
            arr[part[neighbor]] += 1
        arr.append(destNodesCount - sum(arr))
        partEntropy += entropy(arr)
    return partEntropy / len(nodesInPartS)

def findPartitionToSplit(g, nodesInPartS, nodesInPartD, part, sourceNodes):
    maxEntropy = 0
    maxPartition = 0
    for partIndex in range(len(nodesInPartS)):
        avgEntropy = averageEntropy(g, nodesInPartS, nodesInPartD, part, partIndex, sourceNodes)
#         print "average entropy for partition", partIndex, "=", avgEntropy
        if avgEntropy > maxEntropy:
            maxEntropy = avgEntropy;
            maxPartition = partIndex
    return maxPartition

def reGroup(g, nodesInPartS, nodesInPartD, part, sourceNodes):
    sourceNodesCount = len(sourceNodes)
    destNodesCount = len(g.vs) - len(sourceNodes)
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
    

def searchKL(g, nodesInPartS, nodesInPartD, part, sourceNodes, verbose = False):
    changed = True
    while (changed):
        changed = False

        encodingCostDecreased = True
        while(encodingCostDecreased):
            encodingCostDecreased = False
            partIndex = findPartitionToSplit(g, nodesInPartS, nodesInPartD, part, sourceNodes)
            if len(nodesInPartS[partIndex]) > 1:
                currentAverageEntropy = averageEntropy(g, nodesInPartS, nodesInPartD, part, partIndex, sourceNodes)
#                 print "current average entropy", currentAverageEntropy
                for s in nodesInPartS[partIndex]:
                    newNodesInPartS = nodesInPartS[:]
                    newNodesInPartS[partIndex].remove(s)
                    newNodesInPartS.append([s])
                    newPart = part[:]
                    newPart[s] = len(nodesInPartS)
                    newAverageEntropy = averageEntropy(g, newNodesInPartS, nodesInPartD, newPart, partIndex, sourceNodes)
#                     print "new average entropy", newAverageEntropy
                    if newAverageEntropy < currentAverageEntropy:
                        nodesInPartS = newNodesInPartS
                        currentAverageEntropy = newAverageEntropy
                        part = newPart
                if (verbose):
                    print "after split:", nodesInPartS
                reGroup(g, nodesInPartS, nodesInPartD, part, sourceNodes)
                if (verbose):
                    print "after update:", nodesInPartS
                
        merged = True
        while(merged):
            merged = False
            currentTotalCost = totalCost(g, nodesInPartS, nodesInPartD, part, sourceNodes)
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
                    newTotalCost = totalCost(g, newNodesInPartS, nodesInPartD, newPart, sourceNodes)
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

def createPartitionToPartition(g, nodesInPartS, nodesInPartD, part, sourceNodes):
    partitionToPartition = np.zeros(shape=(len(nodesInPartS), len(nodesInPartD)))
    for nodeIndex in range(len(sourceNodes)):
        for neighbor in g.neighbors(sourceNodes[nodeIndex]):
            partitionToPartition[part[sourceNodes[nodeIndex]]][part[neighbor]] += 1
    destNodesCount = len(g.vs) - len(sourceNodes)
    possibleEdgesFromPartitions = np.array([len(elem) for elem in nodesInPartS]) * destNodesCount
    existingEdgesFromPartitions = np.apply_along_axis( sum, axis=1, arr=partitionToPartition)
    complementaryColumnForPartitions = possibleEdgesFromPartitions - existingEdgesFromPartitions
    partitionToPartition = np.append(partitionToPartition, np.transpose(np.matrix(complementaryColumnForPartitions)), axis = 1)
    return partitionToPartition

def initializeDestPartition(destNodes, part):
    destNodePartition = 0
    nodesInPartD = []
    for destNode in destNodes:
        nodesInPartD.append([destNode])
        part[destNode] = destNodePartition
        destNodePartition += 1
    return [part, nodesInPartD]

def partitionGraph(g, iterations):
    totalNodes = len(g.vs)
    part = [0] * totalNodes
    nodesInPartS = []
    sourceNodes = [i for i, x in enumerate(g.vs["type"]) if x == False]
    destNodes = [i for i, x in enumerate(g.vs["type"]) if x == True]
    nodesInPartS.append(sourceNodes[:])
     
    #initialize destPartition so that each node is in its own partition
    initialized = initializeDestPartition(destNodes, part)
    part = initialized[0]
    nodesInPartD = initialized[1]
    
    for i in range(iterations):
        # searchKL for source nodes
        result = searchKL(g, nodesInPartS, nodesInPartD, part, sourceNodes)
        nodesInPartS = result[0]
        part = result[1]
        if (i == 0):
            nodesInPartD = [destNodes[:]]
            for destNode in destNodes:
                part[destNode] = 0
      
        # searchKL for dest nodes    
        result = searchKL(g, nodesInPartD, nodesInPartS, part, destNodes)
        nodesInPartD = result[0]
        part = result[1]
        print "iteration", i, ":\n source:", nodesInPartS, "\n dest:", nodesInPartD, "\n partitioning", part
    return [nodesInPartS, nodesInPartD, part]

# def segmentGraphToPartitions(segment, newGraph):
#     
    

