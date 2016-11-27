package main.java;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by christinebaek on 11/26/16.
 */
public abstract class Align {

    protected final String referencePath;
    protected final String othersPath;
    protected final Boolean horInit;
    protected final Boolean verInit;
    protected final Boolean horEnd;
    protected final Boolean verEnd;
    protected Sequence ref;
    protected List<Sequence> others = new ArrayList<>();
    protected Cell[][] cells;
    protected List<Alignment> alignments = new ArrayList<>();
    protected ScoringMatrix matrix;

    protected Boolean verGapOpen = false;
    protected Boolean horGapOpen = false;
    protected final int gapPenalty = -1;
    protected final int affineGapPenalty = -10;

    public Align(ScoringMatrix matrix, String othersPath, String referencePath, Boolean horInit, Boolean verInit, Boolean horEnd, Boolean verEnd) {
        System.out.println("Alignment Start");
        this.matrix = matrix;
        this.othersPath = othersPath;
        this.referencePath = referencePath;
        this.horInit = horInit;
        this.verInit = verInit;
        this.horEnd = horEnd;
        this.verEnd = verEnd;
        readRef();
        readOthers();
        others.forEach(seq -> {
            System.out.println("*Sequence being processed : " + seq.name);
            align(ref, seq);
        });
        //alignments.forEach(aln -> System.out.println(aln));
        System.out.println("Alignment Analysis Done");
    }

    List<Alignment> getAlignments() {
        return alignments;
    }

    Sequence getRefSeq() {
        return ref;
    }

    /**
     * reads raw fa sequence and parses it
     * @param raw string of the input
     * @return parsed Sequence
     */
    protected Sequence readFa(List<String> raw) {

        Sequence parsed = new Sequence(raw);
        return parsed;
    }

    /**
     * read from reference sequence
     */
    protected void readRef() {
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
    protected void readOthers() {
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

    protected void align(Sequence reference, Sequence other) {
        System.out.println(this.getLocal());
        int refLength = reference.size + 1;
        int otherLength = other.size + 1;
        initMatrix(refLength, otherLength);
        for (int i = 1; i < refLength; i++) {
            for (int j = 1; j < otherLength; j++) {
                char verLetter = reference.sequence.charAt(i-1);
                char horLetter = other.sequence.charAt(j-1);
                getPreviousCell(i, j, verLetter, horLetter, getLocal());
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
    public void getPreviousCell(int i, int j, char verLetter, char horLetter, Boolean local) {
        int verGap = verGapOpen ? gapPenalty : affineGapPenalty;
        int horGap = horGapOpen ? gapPenalty : affineGapPenalty;
        int verticalGapScore = cells[i - 1][j].score + verGap;
        int horizontalGapScore = cells[i][j - 1].score + horGap;

        // this section for diagonal
        int verticalResidueIndex = matrix.residues.indexOf(verLetter);
        int horizontalResidueIndex = matrix.residues.indexOf(horLetter);
        int matchScore = matrix.matrix[verticalResidueIndex][horizontalResidueIndex];
        int diagonalScore = cells[i - 1][j - 1].score + matchScore;

        if (horizontalGapScore < 0 && verticalGapScore < 0 && diagonalScore < 0 && local) {
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
     * initialize matrix, according to type of alignment
     * @param refLength reference sequence length
     * @param otherLength other sequence length
     */
    public abstract void initMatrix(int refLength, int otherLength);

//    /**
//     * update current cell based on DP
//     * @param i row
//     * @param j col
//     * @param verLetter residue of vertical sequence
//     * @param horLetter residue of horizontal sequence
//     */
//    public abstract void getPreviousCell(int i, int j, char verLetter, char horLetter);

    /**
     * find the starting point of traceback
     * @param refLength length of vertical sequence
     * @param otherLength length of horizontal sequence
     * @return
     */
    public abstract Cell findTraceBack(int refLength, int otherLength);

    /**
     *
     * @return whether this is a local alignment
     */
    public abstract Boolean getLocal();
}
