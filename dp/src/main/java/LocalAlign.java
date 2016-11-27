package main.java;

/**
 * Created by christinebaek on 11/26/16.
 */
public class LocalAlign extends Align {

    public LocalAlign(ScoringMatrix matrix, String othersPath, String referencePath, Boolean horInit, Boolean verInit, Boolean horEnd, Boolean verEnd) {
        super(matrix, othersPath, referencePath, horInit, verInit, horEnd, verEnd);
    }

    public void initMatrix(int refLength, int otherLength) {
        cells = new Cell[refLength][otherLength];
        for (int i = 0; i < refLength; i++) {
            for (int j = 0; j < otherLength ; j++) {
                cells[i][j] = new Cell(i, j);
                cells[i][0].score = 0;
                cells[0][j].score = 0;
            }
        }
    }

    /**
     * update current cell based on DP
     * @param i row
     * @param j col
     * @param verLetter residue of vertical sequence
     * @param horLetter residue of horizontal sequence
     */
    public void getPreviousCell(int i, int j, char verLetter, char horLetter) {
        int verGap = verGapOpen ? gapPenalty : affineGapPenalty;
        int horGap = horGapOpen ? gapPenalty : affineGapPenalty;
        int verticalGapScore = cells[i - 1][j].score + verGap;
        int horizontalGapScore = cells[i][j - 1].score + horGap;

        // this section for diagonal
        int verticalResidueIndex = matrix.residues.indexOf(verLetter);
        int horizontalResidueIndex = matrix.residues.indexOf(horLetter);
        int matchScore = matrix.matrix[verticalResidueIndex][horizontalResidueIndex];
        int diagonalScore = cells[i - 1][j - 1].score + matchScore;

        if (horizontalGapScore < 0 && verticalGapScore < 0 && diagonalScore < 0) {
            cells[i][j].score = 0;
            horGapOpen = false;
            verGapOpen = false;
        } else if (horizontalGapScore > verticalGapScore && horizontalGapScore > diagonalScore) {
            cells[i][j].prev = cells[i][j-1];
            cells[i][j].score = horizontalGapScore;
            horGapOpen = true;
            verGapOpen = false;
        } else if (verticalGapScore > diagonalScore && verticalGapScore > horizontalGapScore) {
            cells[i][j].prev = cells[i-1][j];
            cells[i][j].score = verticalGapScore;
            verGapOpen = true;
            horGapOpen = false;
        } else { //diagonalScore is biggest, or there is tie
            cells[i][j].prev = cells[i-1][j-1];
            cells[i][j].score = diagonalScore;
            horGapOpen = false;
            verGapOpen = false;
        }

    }

    /**
     * find the starting point of traceback
     * @param refLength length of vertical sequence
     * @param otherLength length of horizontal sequence
     * @return
     */
    public Cell findTraceBack(int refLength, int otherLength) {
        Cell tbPoint = cells[0][0];
        int bestScore = tbPoint.score;
        for (int i = 0; i < refLength; i++) {
            for (int j = 0; j < otherLength; j++) {
                int curScore = cells[i][j].score;
                if (curScore > bestScore) {
                    tbPoint = cells[i][j];
                    bestScore = curScore;
                }
            }

        }
        return tbPoint;
    }
}

