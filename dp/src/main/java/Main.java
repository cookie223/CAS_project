package main.java;

import java.util.List;

public class Main {

    public static void main(String[] args) throws Exception {
        if (args.length != 4) {
            throw new IllegalArgumentException("incorrect number of arguments");
        }

        boolean horInit = Boolean.parseBoolean(args[0]);
        boolean verInit = Boolean.parseBoolean(args[1]);
        boolean horEnd = Boolean.parseBoolean(args[2]);
        boolean verEnd = Boolean.parseBoolean(args[3]);

        String otherSeq = "src/main/resources/complete_amino_acids.fa";
        //String others = "src/main/resources/file2.fa";
        String referenceSeq = "src/main/resources/s_pyogenes_cas9.fa";


        String matrixDirectory = "src/main/resources/BLOSUM/BLOSUM62";
        ScoringMatrix scoringMatrix = new ScoringMatrix(matrixDirectory);

        String semiAlignFileName = "Cas_SemiGlobal_Alignment";
        SemiGlobalAlign align = new SemiGlobalAlign(scoringMatrix, otherSeq, referenceSeq, horInit, verInit, horEnd, verEnd);
        Sequence refSeq = align.getRefSeq();
        List<Alignment> semiAlign = align.getAlignments();
        new AlignmentToCSV(refSeq, semiAlignFileName, semiAlign);


        String localAlignFileName = "Cas_Local_Alignment";
        LocalAlign local = new LocalAlign(scoringMatrix, otherSeq, referenceSeq, horInit, verInit, horEnd, verEnd);
        Sequence refSequ = local.getRefSeq();
        List<Alignment> localAln = local.getAlignments();
        new AlignmentToCSV(refSequ, localAlignFileName, localAln);


        // test case - comment out for actual run
//        String test = "src/main/resources/test.fa";
//        SemiGlobalAlign testSM = new SemiGlobalAlign(scoringMatrix, test, referenceSeq, horInit, verInit, horEnd, verEnd);
//        Sequence refSeq = testSM.getRefSeq();
//        List<Alignment> semiAlign = testSM.getAlignments();
//        new AlignmentToCSV(refSeq, "testing", semiAlign);



    }
}

