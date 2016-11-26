package main.java;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * aligns many sequences in 'others' to a single reference sequence
 * currently ignores any ties, and prefer diagonal matching
 * Created by christinebaek on 11/25/16.
 */
public class SemiGlobal {

    private final String referencePath;
    private final String othersPath;
    private final Boolean horInit;
    private final Boolean verInit;
    private final Boolean horEnd;
    private final Boolean verEnd;
    private Sequence ref;
    private List<Sequence> others = new ArrayList<>();
    private Cell[][] cells;
    private List<Alignment> alignments = new ArrayList<>();
    private ScoringMatrix matrix;
    private final int endGapPenalty = -3;


    public SemiGlobal(ScoringMatrix matrix, String othersPath, String referencePath, Boolean horInit, Boolean verInit, Boolean horEnd, Boolean verEnd) {
        this.matrix = matrix;
        this.othersPath = othersPath;
        this.referencePath = referencePath;
        this.horInit = horInit;
        this.verInit = verInit;
        this.horEnd = horEnd;
        this.verEnd = verEnd;
        readRef();
        readOthers();
        others.forEach(seq -> align(ref, seq));
        alignments.forEach(aln -> System.out.println(aln));
    }

    /**
     * reads raw fa sequence and parses it
     * @param raw string of the input
     * @return parsed Sequence
     */
    private Sequence readFa(List<String> raw) {

        Sequence parsed = new Sequence(raw);
        return parsed;
    }

    /**
     * read from reference sequence
     */
    private void readRef() {
        String newLine;
        List<String> read = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(referencePath))) {
            while ((newLine = br.readLine()) != null) {read.add(newLine);}
        } catch (FileNotFoundException e) {
            System.out.println("File doesn't exist : " + referencePath);
        } catch (IOException e) {
            System.out.println("Could not read file : " + e);
        }
        ref = readFa(read);
    }

    /**
     * read from non-reference sequence
     */
    private void readOthers() {
        String newLine;
        List<String> read = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(othersPath))) {
            while ((newLine = br.readLine()) != null) {
                if (newLine.contains(">")) {
                    if (read.size() > 0) {
                        others.add(readFa(read));
                    }
                    read = new ArrayList<>();
                }
                read.add(newLine);
            }
        } catch (FileNotFoundException e) {
            System.out.println("File doesn't exist : " + othersPath);
        } catch (IOException e) {
            System.out.println("Could not read file : " + e);
        }
        others.add(readFa(read));
    }

    private void initMatrix(int refLength, int otherLength) {
        cells = new Cell[refLength][otherLength];
        for (int i = 0; i < refLength; i++) {
            for (int j = 0; j < otherLength ; j++) {
                cells[i][j] = new Cell(i, j);
                if (horInit) {
                    cells[i][0].score = endGapPenalty * i;
                    cells[0][j].score = 0;
                } else if (verInit) {
                    cells[i][0].score = 0;
                    cells[0][j].score = endGapPenalty * i;
                }
            }
        }
    }

    private void align(Sequence reference, Sequence other) {
        int refLength = reference.size + 1;
        int otherLength = other.size + 1;
        initMatrix(refLength, otherLength);
        for (int i = 1; i < refLength; i++) {
            for (int j = 1; j < otherLength; j++) {
                char verLetter = reference.sequence.charAt(i-1);
                char horLetter = other.sequence.charAt(j-1);
                getPreviousCell(i, j, verLetter, horLetter);
            }
        }
        Cell tbPoint = findTraceBack(refLength, otherLength);
        Alignment newAln = new Alignment(reference, other, cells, tbPoint);
        alignments.add(newAln);
    }

    /**
     * update current cell based on DP
     * @param i row
     * @param j col
     * @param verLetter residue of vertical sequence
     * @param horLetter residue of horizontal sequence
     */
    private void getPreviousCell(int i, int j, char verLetter, char horLetter) {
        int verticalGapScore = endGapPenalty + cells[i-1][j].score;
        int horizontalGapScore = endGapPenalty + cells[i][j-1].score;

        // this section for diagonal
        int verticalResidueIndex = matrix.residues.indexOf(verLetter);
        int horizontalResidueIndex = matrix.residues.indexOf(horLetter);
        int matchScore = matrix.matrix[verticalResidueIndex][horizontalResidueIndex];
        int diagonalScore = cells[i-1][j-1].score + matchScore;

        if (horizontalGapScore > verticalGapScore && horizontalGapScore > diagonalScore) {
            cells[i][j].prev = cells[i][j-1];
            cells[i][j].score = horizontalGapScore;

        } else if (verticalGapScore > diagonalScore && verticalGapScore > horizontalGapScore) {
            cells[i][j].prev = cells[i-1][j];
            cells[i][j].score = verticalGapScore;
        } else { //diagonalScore is biggest, or there is tie
            cells[i][j].prev = cells[i-1][j-1];
            cells[i][j].score = diagonalScore;
        }

    }

    /**
     * find the starting point of traceback
     * @param refLength length of vertical sequence
     * @param otherLength length of horizontal sequence
     * @return
     */
    private Cell findTraceBack(int refLength, int otherLength) {
        Cell tbPoint = cells[refLength-1][otherLength-1];
        int bestScore = tbPoint.score;
        if (horEnd) {
            // find maximum on the last row
            for (int j = 0; j < otherLength; j++) {
                int curScore = cells[refLength-1][j].score;
                if (curScore > bestScore) {
                    tbPoint = cells[refLength-1][j];
                    bestScore = curScore;
                }
            }
        } else if (verEnd) {
            // find maximum on the last col
            for (int i = 0; i < refLength; i++) {
                int curScore = cells[i][otherLength-1].score;
                if (curScore > bestScore) {
                    tbPoint = cells[i][otherLength-1];
                    bestScore = curScore;
                }
            }
        }
        return tbPoint;
    }

}
