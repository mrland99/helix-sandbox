import numpy as np
import sys

def neighborJoining(D: np.ndarray, n: int, vertexList = None):
    """
        Finds the simple fitting tree for an n x n distance matrix D.

        Input: D = n x n distance matrix, n = number of leaves in D,
                vertexList = maps index in D to vertex in tree
        Returns: tree = adjacency list of tree, treeWeights = edgeWeights between vertices in our tree,
            vertexList = mapping of index of D to vertex in tree
        """

    if vertexList == None:
        vertexList = list(range(n))

    if n == 2:
        tree = {}
        i = vertexList[0]
        j = vertexList[1]
        tree[i] = [j]
        tree[j] = [i]

        edgeLength = D[0][1]
        treeWeights = {}
        treeWeights[(i, j)] = edgeLength
        treeWeights[(j, i)] = edgeLength
        return tree, treeWeights, vertexList

    # construct neighbor joining matrix according to selection function
    # u_i = sum_k D(i, k)
    u = {}
    for i in range(n):
        sum = 0
        for k in range(n):
            sum += D[i][k]
        u[vertexList[i]] = sum
    minScore = sys.maxsize
    i_min = -1
    j_min = -1
    for i in range(n):
        for j in range(i+1, n):
            score = (n-2)*D[i][j]-u[vertexList[i]]-u[vertexList[j]]
            if score < minScore:
                minScore = score
                i_min = i
                j_min = j
    limb_i = 0.5*(D[i_min][j_min] + (u[vertexList[i_min]]-u[vertexList[j_min]])/(n-2))
    limb_j = 0.5*(D[i_min][j_min] + (u[vertexList[j_min]]-u[vertexList[i_min]])/(n-2))

    # add new row column D
    m = vertexList[-1] + 1
    col = np.zeros([n])
    for k in range(n):
        col[k] = 0.5*(D[k][i_min] + D[k][j_min] - D[i_min][j_min])
    row = np.append(col, 0)
    D = np.column_stack((D, col))
    D = np.row_stack((D, row))
    vertexList.append(m)

    # remove rows/columns i and j from D
    D = np.delete(D, (i_min, j_min), 0)
    D = np.delete(D, (i_min, j_min), 1)
    vertex_i = vertexList[i_min]
    vertex_j = vertexList[j_min]
    vertexList.remove(vertex_i)
    vertexList.remove(vertex_j)
    tree, treeWeights, vertexList = neighborJoining(D, n-1, vertexList)

    # Add two new limbs (connecting node m with leaves i and j) to the tree
    if vertex_i not in tree.keys():
        tree[vertex_i] = [m]
    else:
        tree[vertex_i].append(m)
    if vertex_j not in tree.keys():
        tree[vertex_j] = [m]
    else:
        tree[vertex_j].append(m)
    if m not in tree.keys():
        tree[m] = [vertex_i, vertex_j]
    else:
        tree[m].append(vertex_i)
        tree[m].append(vertex_j)

    treeWeights[(vertex_i, m)] = limb_i
    treeWeights[(m, vertex_i)] = limb_i
    treeWeights[(vertex_j, m)] = limb_j
    treeWeights[(m, vertex_j)] = limb_j

    return tree, treeWeights, vertexList