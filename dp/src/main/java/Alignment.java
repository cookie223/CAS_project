package main.java;

/**
 * Created by christinebaek on 11/25/16.
 */
public class Alignment {

    final Sequence reference;
    final Sequence other;
    final int score;
    private final Cell[][] cells;
    String referenceAligned;
    String otherAligned;
    private Cell tracebackPoint;
    private final char gap = "-".charAt(0);
    int endRow;
    int endCol;
    int startRow;
    int startCol;



    public Alignment(Sequence reference, Sequence other, Cell[][] cells, Cell tbp) {
        this.reference = reference;
        this.other = other;
        this.score = tbp.score;
        this.cells = cells;
        this.tracebackPoint = tbp;
        readTraceBack();
    }

    public String toString() {
        return String.format("\nSequence Aligned : %s # with # %s \n" +
                        "Alignment Score : %s \n" +
                        "Reference Subsequence : %s \n" +
                        "Aligned Subsequence : %s \n" +
                        "traceBack started from : [%d][%d] \n" +
                        "traceBack ended at : [%d][%d] \n\n" +
                        "****************************************" +
                        "****************************************"
                , reference.name, other.name, score, referenceAligned, otherAligned, endRow, endCol, startRow, startCol);
    }

    /**
     * trace back to extract the aligned subsequence of the two sequences
     */
    private void readTraceBack() {
        StringBuilder refSB = new StringBuilder();
        StringBuilder othSB = new StringBuilder();

        endRow = tracebackPoint.row;
        endCol = tracebackPoint.col;

        Cell prevCell = new Cell(tracebackPoint.row + 1, tracebackPoint.col + 1);
        Cell curCell = tracebackPoint;
        while (curCell != null) {
            int row = curCell.row - 1;
            int col = curCell.col - 1;
            if (row < 0 || col < 0) break;

            char ver = prevCell.row != curCell.row? reference.sequence.charAt(row) : gap;
            char hor = prevCell.col != curCell.col? other.sequence.charAt(col) : gap;

            refSB.insert(0, ver);
            othSB.insert(0, hor);

            prevCell = curCell;
            curCell = curCell.prev;
        }
        startRow = prevCell.row;
        startCol = prevCell.col;

        referenceAligned = refSB.toString();
        otherAligned = othSB.toString();
    }
}

