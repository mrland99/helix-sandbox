import numpy as np
from sys import maxsize


def fittingAlignment(s1: str, s2: str, scoreMatrix: dict, indelPenalty: int) -> {int, int, str, str}:
    """ Find the highest-scoring fitting alignment of s2 into s1. Note --> s1 > s2
        Input: Two amino acid strings s1, s2. Scoring matrix and indel penalty.
        Return: The maximum alignment score of fitting s2 into s1 followed by an alignment achieving this maximum score.
    """
    # Error checking inputs
    if len(s2) > len(s1):
        print("Error: s1 longer than s2. Please swap the input strings.")
        return -1, -1, -1, -1

    # Initialize Alignment score matrix
    rows = len(s1) + 1
    cols = len(s2) + 1
    matrix = np.zeros([rows, cols], dtype=int)
    for j in range(cols):
        matrix[0][j] = j * indelPenalty

    # Compute Fitting Alignment Matrix
    for i in range(1, rows):
        for j in range(1, cols):
            # check insertion score
            score1 = matrix[i - 1][j] + indelPenalty
            # check match/mismatch score
            key = (s1[i - 1], s2[j - 1])
            if key not in scoreMatrix.keys():
                key = (s2[j - 1], s1[i - 1])
            score2 = matrix[i - 1][j - 1] + scoreMatrix[key]
            # check deletion score
            score3 = matrix[i][j - 1] + indelPenalty

            matrix[i][j] = max(max(score1, score2), score3)

    # Find max alignment score
    maxScore = -maxsize
    maxScoreRow = -1
    maxScoreCol = -1
    for i in range(1, rows):
        if maxScore < matrix[i][cols - 1]:
            maxScore = matrix[i][cols - 1]
            maxScoreRow = i
            maxScoreCol = cols - 1

    # Count the number of optimal alignments
    numOptimal = 0
    for i in range(1, rows):
        if maxScore == matrix[i][cols - 1]:
            numOptimal += 1

    # Find alignment
    v = ""
    w = ""
    i = maxScoreRow
    j = maxScoreCol
    while (i > 0) and (j > 0):
        curr = matrix[i][j]
        left = matrix[i][j - 1]
        top = matrix[i - 1][j]
        diag = matrix[i - 1][j - 1]
        key = (s1[i - 1], s2[j - 1])
        if key not in scoreMatrix.keys():
            key = (s2[j - 1], s1[i - 1])
        if curr == (left + indelPenalty):
            w = s2[j - 1] + w
            v = '-' + v
            j = j - 1
        elif curr == (top + indelPenalty):
            w = '-' + w
            v = s1[i - 1] + v
            i = i - 1
        elif curr == (diag + scoreMatrix[key]):
            w = s2[j - 1] + w
            v = s1[i - 1] + v
            i = i - 1
            j = j - 1
        else:
            print("Error in generating sequence alignment")
            print("i = " + str(i) + ", j = " + str(j))
            break

    return maxScore, numOptimal, v, w
