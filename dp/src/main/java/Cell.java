package main.java;

/**
 * Created by christinebaek on 11/25/16.
 */
public class Cell {
    int score;
    Cell prev;
    int row;
    int col;

    public Cell(int row, int col) {
        this.row = row;
        this.col = col;
    }

}
