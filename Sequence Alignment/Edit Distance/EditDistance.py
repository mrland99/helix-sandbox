import numpy as np


def editDistance(s1: str, s2: str) -> {int}:
    """ The edit distance is the minimum number of edit operations needed to transform s1 into s2, where operation is
        defined as substitution, insertion, or deletion of a single symbol.
        Input: Two amino acid strings s1, s2.
        Return: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
    """

    # Initialize score matrix
    rows = len(s1) + 1
    cols = len(s2) + 1
    matrix = np.zeros([rows, cols], dtype=int)
    penalty = 1
    for i in range(rows):
        matrix[i][0] = i * penalty
    for j in range(cols):
        matrix[0][j] = j * penalty

    # Compute Edit Distance Matrix
    for i in range(1, rows):
        for j in range(1, cols):
            # check insertion score
            score1 = matrix[i - 1][j] + penalty
            # check match/mismatch score
            if s1[i - 1] == s2[j - 1]:
                score2 = matrix[i-1][j-1]
            else:
                score2 = matrix[i-1][j-1] + 1
            # check deletion score
            score3 = matrix[i][j - 1] + penalty

            matrix[i][j] = min(min(score1, score2), score3)

    dist = matrix[rows - 1][cols - 1]
    return dist
