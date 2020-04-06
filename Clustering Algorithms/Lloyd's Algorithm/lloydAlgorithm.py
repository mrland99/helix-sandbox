import numpy as np
import sys
import math


def LloydAlgorithm(centroids: np.ndarray, points: np.ndarray):
    """
    Lloyd's algorithm is a popular clustering heuristic for k-Means clustering problem.
    It first chooses k arbitrary points Centers from Data as centers and then iteratively performs the following
     two steps:

        1) Centers to Clusters: After centers have been selected, assign each data point to the cluster corresponding
        to its nearest center; ties are broken arbitrarily.
        2) Clusters to Centers: After data points have been assigned to clusters, assign each cluster’s center of
        gravity to be the cluster’s new center.

    :param centroids: initial centers
    :param points: numpy array of data points
    :return: set of centers consisting of k points
    """
    # store number of centers and dimension
    k, m = centroids.shape

    labels = clusterPoints(centroids, points)
    newCentroids = updateCentroids(points, labels, k, m)
    while not np.allclose(centroids, newCentroids, atol=1e-03):
        centroids = newCentroids
        labels = clusterPoints(centroids, points)
        newCentroids = updateCentroids(points, labels, k, m)
    return centroids


def clusterPoints(centroids, points):
    """
    assign data points to clusters (labels)
    :param centroids: numpy array of centroids
    :param points: numpy array of points
    :return: list of labels
    """
    labels = []
    for p in range(len(points)):
        dist = sys.maxsize
        label = -1
        for c in range(len(centroids)):
            newDist = distance(centroids[c], points[p])
            if newDist < dist:
                dist = newDist
                label = c
        labels.append(label)
    return labels


def updateCentroids(points, labels, k, m):
    """
    assign each cluster's center of gravity to be the cluster's new center
    :param points: numpy array of data points
    :param labels: list of labels that map each data point to a cluster
    :param k: number of centers/centroids
    :param m: number of dimensions of data
    :return: updated numpy array of centers/centroids
    """
    centroids = np.zeros((k, m))
    for c in range(len(centroids)):
        s = np.zeros((1, m))
        count = 0
        for p in range(len(labels)):
            if labels[p] == c:
                s = np.add(s, points[p])
                count += 1
        centroids[c] = np.true_divide(s, count)
    return centroids


def distance(a, b):
    if len(a) != len(b):
        print("Error: a and b have different dimensions!")
        return -1
    dist = 0
    for i in range(len(a)):
        dist += (a[i] - b[i]) ** 2
    return math.sqrt(dist)
