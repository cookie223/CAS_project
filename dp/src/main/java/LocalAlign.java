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

    public Boolean getLocal() { return true; }

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

