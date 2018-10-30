# gene-sequencing
A dynamic programming algorithm for computing the minimal cost of aligning gene sequences (taxa) and for extracting optimal alignments

In light of the SARS outbreak in 2003, I have chosen to use the SARS virus as one of my DNA sequences. SARS is a coronavirus, and I have also provided the genomes for several other coronaviruses to be compared to the SARS virus. 

### Design
The viruses are compared by calculating their 'edit distances'. Edit distances are weighted using the Needleman/Wunsch cost function.

### Needleman/Wunsch Costs
There are three types of edits when comparing two gene sequences and each has a specific weight. Weights were determined based on what would statistically occur more often in biological changes between two gene sequences. The edit types and their weights are explained  in the table below.

| **Type** | **Weight** | **Description** | **Example** (Indexes start at 0) | 
|----------|------------|-----------------|----------------------------------|
| Match | -3 | When there is no change at a specific index, we like these | Comparing index 1 of 'c**t**g' with 'g**t**c' we see there is a match of 't'. Edit weight of **-3**
| Substitution | 1 | Single character mismatches that could be changed by changing one of the sub-sequences indices | Comparing index 1 of 'c**t**g' with 'c**g**c' we could substitute the 't' in 'ctg' with a 'g'. Edit weight of **1**
| Insert/Delete (indel) | 5 | These are bad (and don't happen often in biology. They essentially create gaps due to there not being enough characters in one sequence | Comparing index 3 of 'ctg' with 'ctg**t**c' we would need to insert 't' and 'c' to the end of 'ctg'. Edit weight of **10**

### Implementation
There are two algorithms for finding the optimal alignments:
1. The first involves creating an NxN matrix where N is the specified length of each gene to compare. Then comparing both sub-sequences and evaluating their differences. This method will provide the most optimal answer but is also a O(n^2) algorithm. The cell in the matrix at position N-1, N-1 (the bottom right cell) will give us the optimal 'edit distance' between the two gene sequences.
2. The second is a 'banded' method. This algorithm esstentially only worries about the diagonal of the same NxN matrix from the previous algorithm given a limit to how many indels can occur in a row before we decide to exclude the comparisons beyond that point. The 'band' is given a width that represents the limit of indels.  If the bottom right cell, mentioned above, is not present in the banded matrix then we discard the comparison of the two gene sequences because they are so different that we don't care. We do lose some optimization with this algorithm but it is a O(n) algorithm.
