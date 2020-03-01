# helix-sandbox
A compilation of various useful algorithms in computational biology written by myself. Work in progress.

Here are the descriptions of the algorithms:

1. Sequence Alignment:
   * `Edit Distance` - Given two strings, compute the minimum number of insertions, deletions, and substitutions required to transform one string into the other.
   * `Global Alignment`- Aligns two sequences from beginning to end factoring in mistmatches and indels. Each letter in each sequence is aligned only once.
   * `Local Alignment` - Given sequence *v* and *w*, returns optimal alignment between all substrings of *v* and *w*. Another way to think about it is maximum global alignment score over all substring pairs of *v* and *w*.
   * `Overlap Alignment` - Also known as `Semi-Global Alignment`. It is alignment where start and end gaps (gaps appearing as prefixes or suffixes) are ignored.
   * `Fitting Alignment` - Given sequence *v* and *w*, returns optimal alignment between the whole sequence *w* to all possible substrings of *v*. 
2. Phylogeny:
   * `Additive Phylogeny` - Finds the simple tree fitting an n x n additive distance matrix D.
   
