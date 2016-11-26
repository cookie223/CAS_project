package main.java;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * parses NCBI standard format PAM/BLOSUM matrix
 * Created by christinebaek on 11/25/16.
 */
public class ScoringMatrix {

    private static final int SIZE = 24;
    private final String path;
    String name; // name of matrix
    List<Character> residues; // name of residues
    int[][] matrix; // scoring matrix


    public ScoringMatrix(String path) {

        this.path = path;
        List<String> raw = readFile(this.path);
        this.name = toString();
        readMatrix(raw);
    }

    @Override
    public String toString() {
        int index = path.lastIndexOf("/") + 1;
        return path.substring(index);
    }

    /**
     * read from scoring matrix file
     * @param path of the file
     * @return raw string of the scores
     */
    private List<String> readFile(String path) {
        String newLine;
        List<String> read = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            while ((newLine = br.readLine()) != null) {
                if (!newLine.contains("#") && !newLine.startsWith(" "))
                read.add(newLine);}
        } catch (FileNotFoundException e) {
            System.out.println("File doesn't exist : " + path);
        } catch (IOException e) {
            System.out.println("Could not read file : " + e);
        }

        return read;
    }

    private void readMatrix(List<String> raw) {
        matrix = new int[SIZE][SIZE];
        residues = new ArrayList<>();
        for (int i = 0; i < raw.size(); i++) {
            String[] split = raw.get(i).split("\\s+");
            for (int j = 0; j < split.length; j++) {
                if (j == 0) {
                    residues.add(split[j].charAt(0));
                } else {
                    matrix[i][j-1] = Integer.parseInt(split[j]);
                }
            }
        }
    }

}
