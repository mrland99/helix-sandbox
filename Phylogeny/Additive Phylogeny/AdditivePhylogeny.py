import numpy as np
import math

def additivePhylogeny(D: np.ndarray, n: int, newIndexCount: int) -> (dict, dict):
    """
    Finds the simple fitting tree for an n x n distance matrix D.

    Input: D = distance matrix, n = number of leaves in D, newIndexCount = number of original leaves
    Returns: tree = adjacency list of tree, treeWeights = edgeWeights between vertices in our tree
    """
    if n == 2:
        tree = {}
        tree[0] = [1]
        tree[1] = [0]

        edgeLength = D[0][1]
        treeWeights = {}
        treeWeights = {}
        treeWeights[(0, 1)] = edgeLength
        treeWeights[(1, 0)] = edgeLength
        return tree, treeWeights, newIndexCount

    limbLength = getLimbLength(D, n)
    for j in range(1, n):
        D[j-1][n-1] = D[j-1][n-1] - limbLength
        D[n-1][j-1] = D[j-1][n-1]

    leaf_i_index = -1
    leaf_k_index = -1
    x = -1
    stop = False
    for i in range(n):
        for k in range(n):
            if (D[i][k] == D[i][n-1] + D[n-1][k]):
                x = D[i][n-1]
                leaf_i_index = i
                leaf_k_index = k
                stop = True
                break
        if stop:
            break
    x = D[leaf_i_index][n-1]
    tree, treeWeights, newIndexCount = additivePhylogeny(D[:-1, :-1], n-1, newIndexCount)
    # Add new leaves to tree and update accordingly
    nearVertex1, nearVertex2, dist1, dist2 = nearest(tree, treeWeights, x, leaf_i_index, leaf_k_index)
    newVertex = nearVertex1
    # Check to make sure we need to create new node instead of appending it to end
    if dist1 != 0:
        newVertex = newIndexCount
        newIndexCount += 1
        # insert newVertex
        tree[nearVertex1].remove(nearVertex2)
        tree[nearVertex1].append(newVertex)
        treeWeights[(nearVertex1, newVertex)] = dist1
        treeWeights[(newVertex, nearVertex1)] = dist1

        tree[nearVertex2].remove(nearVertex1)
        tree[nearVertex2].append(newVertex)
        treeWeights[(nearVertex2, newVertex)] = dist2
        treeWeights[(newVertex, nearVertex2)] = dist2

        tree[newVertex] = [nearVertex1, nearVertex2]
        del treeWeights[(nearVertex1, nearVertex2)]
        del treeWeights[(nearVertex2, nearVertex1)]

    # add limb
    tree[newVertex].append(n-1)
    tree[n-1] = [newVertex]
    treeWeights[(newVertex, n-1)] = limbLength
    treeWeights[(n-1, newVertex)] = limbLength

    return tree, treeWeights, newIndexCount

def nearest(tree, treeWeights, x, i, k) -> (int, int, int, int):
    """
        Finds the nearest two leaves on the path i -> k to our new vertex
        Return: two nearest vertices to our new vertex and their respective weights
    """
    # BFS to find path between i and k
    queue = [[i]]
    visited = set([i])
    findPath = []
    while len(queue) > 0:
        path = queue.pop()
        vertex = path[-1]
        visited.add(vertex)
        if vertex == k:
            findPath = path
            break
        else:
            for nextVertex in tree[vertex]:
                if nextVertex not in visited:
                    queue.append(path + [nextVertex])
    # Determine which two nodes the new vertex lies between
    dist = 0
    for k in range(len(findPath) - 1):
        i, j = findPath[k], findPath[k+1]
        if (dist + treeWeights[(i, j)] > x):
            return (i, j, x - dist, dist + treeWeights[(i, j)] - x)
        dist += treeWeights[(i, j)]

def getLimbLength(D: np.ndarray, n: int) -> int:
    """ Given an additive distance matrix, returns the limb length of leaf n.
            Input: D = distance matrix, n = leaf we are calculating limb length of
            Return: length of limb.
        """

    # leaf n corespondes to index n - 1
    j = n - 1
    limbLength = math.inf;
    if len(D) < 3:
        print("Error: Not enough leaves in tree (less than 2).")
        return -1

    # Fix arbitrary i
    i = j - 1
    if i < 0:
        i = j + 1

    for k in range(len(D)):
        if k != i and k != j:
            length = (D[i][j] + D[j][k] -D[i][k])/2
            if length < limbLength:
                limbLength = length

    return int(limbLength)