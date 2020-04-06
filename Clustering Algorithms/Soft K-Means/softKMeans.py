import math
import numpy as np


def softCluster(B: float, centroids: np.ndarray, points: np.ndarray, n=100):
    """
    The soft k-means clustering algorithm starts from randomly chosen centers and iterates the following two steps:

        1) Centers to Soft Clusters (E-step): After centers have been selected, assign each data point a
        “responsibility” value for each cluster, where higher values correspond to stronger cluster membership.
        2) Soft Clusters to Centers (M-step): After data points have been assigned to soft clusters, compute
        new centers.

    :param n: number of times e-step and m-step are run
    :param B: stiffness parameter used calculation later on
    :param centroids: initial centers
    :param points: data points
    :return: centers after running the algorithm n steps
    """
    for i in range(n):
        hiddenMatrix = e_step(B, centroids, points)
        centroids = m_step(hiddenMatrix, points)

    return centroids


def m_step(hiddenMatrix, points):
    k = len(hiddenMatrix)
    m = points.shape[1]
    newCentroids = np.zeros((k, m))
    for i in range(k):
        for j in range(m):
            newCentroids[i][j] = calculate_m_step_entry(i, j, hiddenMatrix, points)
    return newCentroids


def calculate_m_step_entry(i, j, hiddenMatrix, points):
    numerator = np.matmul(hiddenMatrix[i], points[:, j])
    denominator = np.matmul(hiddenMatrix[i], np.transpose(np.ones(len(points))))
    return numerator / denominator


def e_step(B, centroids, points):
    k = len(centroids)
    n = len(points)
    hiddenMatrix = np.zeros((k, n))
    for i in range(k):
        for j in range(n):
            hiddenMatrix[i][j] = calculate_e_step_entry(i, j, B, centroids, points)
    return hiddenMatrix


def calculate_e_step_entry(i, j, B, centroids, points):
    val = -1
    numerator = math.exp(-B * distance(centroids[i], points[j]))
    denominator = 0
    for k in range(len(centroids)):
        denominator += math.exp(-B * distance(centroids[k], points[j]))
    return numerator / denominator


def distance(a, b):
    if len(a) != len(b):
        print("Error: a and b have different dimensions!")
        return -1
    dist = 0
    for i in range(len(a)):
        dist += (a[i] - b[i]) ** 2
    return math.sqrt(dist)
