import numpy as np

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

def writePartitionedGraphToFile(nodesInPartS, nodesInPartD, part, adj, sourceNodes, destNodes, fileName):
    writeNumbersToFile([len(sourceNodes), len(destNodes), len(nodesInPartS), len(nodesInPartD)], fileName)
    writeArrayToFile(part, fileName)
    for i in range(len(nodesInPartS)):
        nodesInPartS[i].sort()
        for j in range(len(nodesInPartD)):
            nodesInPartD[j].sort()
            adj2 = [[adj[row][col] for col in nodesInPartD[j]] for row in nodesInPartS[i]]
            writeMatrixToFile(adj2, fileName)
            
def writeInitialGraphToFile(adj, sourceNodes, destNodes, fileName):
    writeNumbersToFile([len(sourceNodes), len(destNodes)], fileName)
    adj2 = [[adj[row][col] for col in destNodes] for row in sourceNodes]
    writeMatrixToFile(adj2, fileName)

def writeMatrixToPnms(adj, fileName1, fileName2, sourceNodes, destNodes, nodesInPartS, nodesInPartD):
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