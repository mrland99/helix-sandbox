import numpy as np

def loadMatrixFile(filename):
    """
    Reads in text as n and a tab-delimited n x n additive matrix.
    Sample File:
        4
        0   13  21  22
        13  0   12  13
        21  12  0   13
        22  13  13  0
    """
    with open(filename) as file:
        dim = int(file.readline())
        matrix = [[int(entry) for entry in line.split()] for line in file]
        matrix = np.array([np.array(row) for row in matrix])
        return matrix, dim

def printPhylogenyTree(tree: dict, treeWeights: dict):
    """
    Prints Phylogeny tree in readable format. The n leaves are labeled 0, 1, ..., n-1 in order of their appearance
    in the distance matrix.

    Input: tree = adjacency list of tree, treeWeights = edgeWeights between vertices in our tree
    """
    for i in sorted(tree):
        for j in sorted(tree[i]):
            print(str(i) + "->" + str(j) + ":" + str(treeWeights[(i,j)]))