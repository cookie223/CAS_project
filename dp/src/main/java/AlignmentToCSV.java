package main.java;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by christinebaek on 11/26/16.
 */
public class AlignmentToCSV {

    private Sequence ref;
    private String name;
    private List<Alignment> alignments;
    private final int leftMenu = 6;
    private static final int COLSIZE = 10;

    public AlignmentToCSV(Sequence ref, String name, List<Alignment> alignments) throws Exception {
        this.ref = ref;
        this.name = name + ".csv";
        this.alignments = alignments;
        toCSV();
    }

    private void toCSV() throws Exception {
        System.out.println(String.format("CSV write starting for %s", name));
        String csvFile = name;
        FileWriter writer = new FileWriter(csvFile);

        // write header
        List<String> positions = new ArrayList<>();

        positions.add("Alignment Score");
        positions.add("Average Score per Base");
        positions.add("Percent Sequence aligned");
        positions.add("Traceback Start");
        positions.add("Traceback End");
        positions.add("Sequence");

        for (int i = leftMenu; i < ref.size/10 + leftMenu; i++) {
            positions.add(String.format("%d", (i-leftMenu) * 10));
        }
        CSVUtils.writeLine(writer, positions);


        // write content
        for(Alignment aln : alignments) {


            List<String> newLine = new ArrayList<>();

            newLine.add(String.format("%s", aln.score));
            newLine.add(String.format("%s", aln.avgScore));
            newLine.add(String.format("%s", aln.percent));
            newLine.add(String.format("[%d][%d]", aln.endRow, aln.endCol));
            newLine.add(String.format("[%d][%d]", aln.startRow, aln.startCol));
            newLine.add(aln.other.name);

            int startPos = aln.startRow/COLSIZE;
            int endPos = aln.endRow/COLSIZE;

            for (int i = leftMenu; i < ref.size/10 + leftMenu; i++) {
                if (i >= startPos && i < endPos) { newLine.add("x"); }
                else { newLine.add(""); }
            }
            CSVUtils.writeLine(writer, newLine);

        }


        // close

        writer.flush();
        writer.close();
    }
}
