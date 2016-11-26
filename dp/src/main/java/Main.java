package main.java;

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
        String referenceSeq = "src/main/resources/s_pyogenes_cas9.fa";


        String matrixDirectory = "src/main/resources/BLOSUM/BLOSUM62";
        ScoringMatrix scoringMatrix = new ScoringMatrix(matrixDirectory);

        SemiGlobal align = new SemiGlobal(scoringMatrix, otherSeq, referenceSeq, horInit, verInit, horEnd, verEnd);


        // test case - comment out for actual run
        String test = "src/main/resources/test.fa";
        //SemiGlobal testSM = new SemiGlobal(scoringMatrix, test, referenceSeq, horInit, verInit, horEnd, verEnd);



    }
}

