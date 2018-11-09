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
                    # If the sub sequences are the same then we have a full match
                    if (subSeq2 == subSeq1):
                        score = len(subSeq2) * MATCH
                        alignment1 = subSeq1.format(
                            i+1, len(sequences[i]), align_length, ',BANDED' if banded else '')
                        alignment2 = subSeq2.format(
                            ',BANDED' if banded else '')
                    else:
                        # Otherwise we need to find the optimal alignment
                        if (self.banded):
                            score = self.alignBanded(subSeq1, subSeq2)
                        else:
                            score = self.alignUnrestricted(subSeq1, subSeq2)
                        # Backtrace for the best edit
                        if (score != math.inf):
                            alignments = self.findOptimalPath(subSeq1, subSeq2)
                            alignment2 = alignments[1].format(
                                ',BANDED' if banded else '')
                            alignment1 = alignments[0].format(
                                i+1, len(sequences[i]), align_length, ',BANDED' if banded else '')
                        else:
                            alignment1 = 'No Alignment Possible'
                            alignment2 = 'No Alignment Possible'

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

        self.initAlignmentArrays(lenSubSeq1, lenSubSeq2)

        for row in range(1, lenSubSeq1 + 1):
            for col in range(1, lenSubSeq2 + 1):
                if (subSeq1[row - 1] == subSeq2[col - 1]):
                    difference = MATCH
                else:
                    difference = SUB
                upperVal = (self.editScores[row - 1, col] + INDEL)
                leftVal = (self.editScores[row, col - 1] + INDEL)
                diagonalVal = (self.editScores[row - 1, col - 1] + difference)
                minVal = min(leftVal, upperVal, diagonalVal)
                self.editScores[row, col] = minVal
                if (minVal == diagonalVal):
                    self.backTrace[row, col] = [DIAGONAL, difference]
                elif (minVal == upperVal):
                    self.backTrace[row, col] = [UP, INDEL]
                else:
                    self.backTrace[row, col] = [LEFT, INDEL]

        return self.editScores[lenSubSeq1, lenSubSeq2]

    def alignBanded(self, subSeq1, subSeq2):
        lenSubSeq1 = len(subSeq1)
        lenSubSeq2 = len(subSeq2)

        lenDiff = abs(lenSubSeq1 - lenSubSeq2)
        # Not possible. Don't even bother trying
        if(lenDiff > MAXINDELS):
            return math.inf

        self.initAlignmentArrays(lenSubSeq1, lenSubSeq2)

        # Offset for fiding the correct index to use while in the middle section of the 'band'
        middleSectStartIdx = MAXINDELS + 1
        middleSectEndIdx = (lenSubSeq1 + 1) - lenDiff - MAXINDELS

        for row in range(1, MAXINDELS + 1):
            for col in range(1, (2 * MAXINDELS) + 1):
                # We are in either the beginning of the 'band', indexing is like unrestricted
                if(row == 1 and col > middleSectStartIdx):
                    break
                if(row == 2 and col > middleSectStartIdx + 1):
                    break

                if (subSeq1[row - 1] == subSeq2[col - 1]):
                    difference = MATCH
                else:
                    difference = SUB
                upperVal = (self.editScores[row - 1, col] + INDEL)
                leftVal = (self.editScores[row, col - 1] + INDEL)
                diagonalVal = (
                    self.editScores[row - 1, col - 1] + difference)
                minVal = min(leftVal, upperVal, diagonalVal)
                self.editScores[row, col] = minVal
                if (minVal == diagonalVal):
                    self.backTrace[row, col] = [DIAGONAL, difference]
                elif (minVal == upperVal):
                    self.backTrace[row, col] = [UP, INDEL]
                else:
                    self.backTrace[row, col] = [LEFT, INDEL]

        # We are in the middle section of the 'band' and the indexing gets a little strange
        strOffset = 0
        for row in range(middleSectStartIdx, middleSectEndIdx + 1):
            for col in range(0, (2 * MAXINDELS) + 1):
                if (subSeq1[row - 1] == subSeq2[col + strOffset]):
                    difference = MATCH
                else:
                    difference = SUB
                diagonalVal = (self.editScores[row - 1, col] + difference)
                upperVal = math.inf
                leftVal = math.inf
                # Avoid index out of bounds
                if(col != 2 * MAXINDELS):
                    upperVal = (self.editScores[row - 1, col + 1] + INDEL)
                # Avoid index out of bounds
                if(col + 1 >= strOffset):
                    leftVal = (self.editScores[row, col - 1] + INDEL)
                minVal = min(leftVal, upperVal, diagonalVal)
                self.editScores[row, col] = minVal
                if (minVal == diagonalVal):
                    self.backTrace[row, col] = [DIAGONAL, difference]
                elif (minVal == upperVal):
                    self.backTrace[row, col] = [UP, INDEL]
                else:
                    self.backTrace[row, col] = [LEFT, INDEL]
            # Increment the index offset to make sure we are comparing the correct subsection of subSeq2
            strOffset += 1

        # for col in range(0, (2 * MAXINDELS) + 1):
        #     if (subSeq1[middleSectEndIdx] == subSeq2[col + strOffset]):
        #         difference = MATCH
        #     else:
        #         difference = SUB
        #     upperVal = (self.editScores[middleSectEndIdx + 1, col] + INDEL)
        #     leftVal = (self.editScores[middleSectEndIdx + 1, col - 1] + INDEL)
        #     diagonalVal = (
        #         self.editScores[middleSectEndIdx + 1, col - 1] + difference)
        #     minVal = min(leftVal, upperVal, diagonalVal)
        #     self.editScores[middleSectEndIdx + 1, col] = minVal
        #     if (minVal == diagonalVal):
        #         self.backTrace[middleSectEndIdx +
        #                        1, col] = [DIAGONAL, difference]
        #     elif (minVal == upperVal):
        #         self.backTrace[middleSectEndIdx + 1, col] = [UP, INDEL]
        #     else:
        #         self.backTrace[middleSectEndIdx + 1, col] = [LEFT, INDEL]

                # We are in the last section
                # Use this column counter as another offset for keeping the infs in the editScores table
        colCounter = 1
        for row in range(middleSectEndIdx + 2, lenSubSeq1 + 1):
            for col in range(colCounter, (2 * MAXINDELS) + 1):
                if (subSeq1[row - 1] == subSeq2[col + strOffset]):
                    difference = MATCH
                else:
                    difference = SUB
                upperVal = (self.editScores[row - 1, col] + INDEL)
                leftVal = (self.editScores[row, col - 1] + INDEL)
                diagonalVal = (
                    self.editScores[row - 1, col - 1] + difference)
                minVal = min(leftVal, upperVal, diagonalVal)
                self.editScores[row, col] = minVal
                if (minVal == diagonalVal):
                    self.backTrace[row, col] = [DIAGONAL, difference]
                elif (minVal == upperVal):
                    self.backTrace[row, col] = [UP, INDEL]
                else:
                    self.backTrace[row, col] = [LEFT, INDEL]
            colCounter += 1

        return self.editScores[lenSubSeq1, 2 * MAXINDELS]

    def initAlignmentArrays(self, lenSubSeq1, lenSubSeq2):
        if(self.banded):
            self.editScores = np.full(
                (lenSubSeq1 + 1, (2 * MAXINDELS) + 1), math.inf)
            self.backTrace = np.ndarray(
                (lenSubSeq1 + 1, (2 * MAXINDELS) + 1), dtype=object)
            self.editScores[0, 0] = 0
            self.backTrace[0, 0] = [START, 0]

            for x in range(1,  MAXINDELS + 1):
                self.editScores[x, 0] = x * INDEL
                self.backTrace[x, 0] = [UP, INDEL]

            for y in range(1,  MAXINDELS + 1):
                self.editScores[0, y] = y * INDEL
                self.backTrace[0, y] = [LEFT, INDEL]
        else:
            self.editScores = np.empty((lenSubSeq1 + 1, lenSubSeq2 + 1))
            self.backTrace = np.ndarray(
                (lenSubSeq1 + 1,  lenSubSeq2 + 1), dtype=object)

            self.editScores[0, 0] = 0
            self.backTrace[0, 0] = [START, 0]

            for x in range(1,  lenSubSeq1 + 1):
                self.editScores[x, 0] = x * INDEL
                self.backTrace[x, 0] = [UP, INDEL]

            for y in range(1, lenSubSeq2 + 1):
                self.editScores[0, y] = y * INDEL
                self.backTrace[0, y] = [LEFT, INDEL]

    def findOptimalPath(self, subSeq1, subSeq2):
        lenSubSeq1 = len(subSeq1)
        lenSubSeq2 = len(subSeq2)
        i = lenSubSeq1
        j = lenSubSeq2
        backwardStr1 = ''
        backwardStr2 = ''
        while(self.backTrace[i, j][0] != START):
            if(self.backTrace[i, j][0] == DIAGONAL):
                i -= 1
                j -= 1
                backwardStr1 += subSeq1[i]
                backwardStr2 += subSeq2[j]

            elif(self.backTrace[i, j][0] == LEFT):
                j -= 1
                backwardStr2 += subSeq2[j]
                backwardStr1 += '-'

            elif(self.backTrace[i, j][0] == UP):
                i -= 1
                backwardStr1 += subSeq1[i]
                backwardStr2 += '-'

        return [backwardStr1[::-1], backwardStr2[::-1]]

    def difference(self, subSeq1, subSeq2, idx1, idx2):
        if (subSeq1[idx1] == subSeq2[idx2]):
            return MATCH
        else:
            return SUB
