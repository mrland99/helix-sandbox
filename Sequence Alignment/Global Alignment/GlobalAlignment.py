import numpy as np


def globalAlignment(s1: str, s2: str, scoreMatrix: dict, indelPenalty: int) -> {int, str, str}:
    """ Find the highest-scoring alignment between two strings using a scoring matrix.
        Input: Two amino acid strings s1, s2. Scoring matrix and indel penalty.
        Return: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
    """

    # Initialize Alignment score matrix
    rows = len(s1) + 1
    cols = len(s2) + 1
    matrix = np.zeros([rows, cols], dtype=int)
    for i in range(rows):
        matrix[i][0] = i * indelPenalty
    for j in range(cols):
        matrix[0][j] = j * indelPenalty

    # Compute Global Alignment Matrix
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

    maxScore = matrix[rows - 1][cols - 1]
    # Find alignment
    v = ""
    w = ""
    i = rows - 1
    j = cols - 1
    while (i > 0) or (j > 0):
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

    return maxScore, v, w
