package main.java;

import java.util.List;

/**
 * Created by christinebaek on 11/25/16.
 */
public class Sequence {

    String name;
    String sequence;
    int size;

    @Override
    public String toString() {
        return sequence;
    }


    Sequence(List<String> raw) {
        parse(raw);
    }

    private void parse(List<String> raw) {
        StringBuilder sb  = new StringBuilder();
        raw.forEach(line->{
            if (line.contains(">")) {
                name = line.substring(1);
            } else {
                sb.append(line);
            }
        });
        sequence = sb.toString();
        size = sequence.length();
    }

}

