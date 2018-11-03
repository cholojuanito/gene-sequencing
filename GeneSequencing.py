#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import numpy as np

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

# Used for backtracing
UP = 'u'
LEFT = 'l'
DIAGONAL = 'd'
START = 's'


class GeneSequencing:

    def __init__(self):
        pass


# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        for i in range(len(sequences)):
            jresults = []
            subSeq1 = sequences[i][:self.MaxCharactersToAlign]
            for j in range(len(sequences)):
                subSeq2 = sequences[j][:self.MaxCharactersToAlign]

                if(j < i):
                    s = {}
                else:
                    ###################################################################################################
                    # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
                    alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                                                                                           len(sequences[i]), align_length, ',BANDED' if banded else '')
                    alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                                                                                           len(sequences[j]), align_length, ',BANDED' if banded else '')

                    if (self.banded):
                        self.alignBanded(subSeq1, subSeq2)
                    else:
                        score = self.alignUnrestricted(subSeq1, subSeq2)

                    s = {'align_cost': score, 'seqi_first100': alignment1,
                         'seqj_first100': alignment2}
                    table.item(i, j).setText('{}'.format(
                        int(score) if score != math.inf else score))
                    table.repaint()
                jresults.append(s)
            results.append(jresults)
        return results

    def alignUnrestricted(self, subSeq1, subSeq2):
        lenSubSeq1 = len(subSeq1)
        lenSubSeq2 = len(subSeq2)
        self.editScores = np.zeros((lenSubSeq1 + 1, lenSubSeq2 + 1))
        self.backTrace = np.ndarray(
            (lenSubSeq1 + 1,  lenSubSeq2 + 1), dtype=object)
        self.backTrace.fill(START)

        for x in range(1,  lenSubSeq1 + 1):
            self.editScores[x, 0] = x * INDEL
            self.backTrace[x, 0] = UP

        for y in range(1, lenSubSeq2 + 1):
            self.editScores[0, y] = y * INDEL
            self.backTrace[0, y] = LEFT

        for row in range(1, lenSubSeq1 + 1):
            for col in range(1, lenSubSeq2 + 1):
                leftVal = (self.editScores[row - 1, col] + INDEL)
                upperVal = (self.editScores[row, col - 1] + INDEL)
                diagonalVal = (self.editScores[row - 1, col - 1] +
                               self.difference(subSeq1, subSeq2, row - 1, col - 1))
                minVal = min(leftVal, upperVal, diagonalVal)
                self.editScores[row, col] = minVal
                if (minVal == diagonalVal):
                    self.backTrace[row, col] = DIAGONAL
                elif (minVal == upperVal):
                    self.backTrace[row, col] = UP
                else:
                    self.backTrace[row, col] = LEFT

        return self.editScores[lenSubSeq1, lenSubSeq2]

    def alignBanded(self, subSeq1, subSeq2):
        pass

    def difference(self, subSeq1, subSeq2, idx1, idx2):
        if (subSeq1[idx1] == subSeq2[idx2]):
            return MATCH
        else:
            return SUB
