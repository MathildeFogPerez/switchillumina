package ch.irb.blast.switchregion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 Copyright 2017 - Mathilde Foglierini Perez
 This code is distributed open source under the terms of the GNU Free Documention License.
 This class is used to parse a blast output file and o create an ArrayList of Alignment objects.
 It works only if we have ONE query sequence (otherwise use SwitchBlastParserPerBlock class)
 */
public class SwitchBlastParser {

    private ArrayList<Alignment> alignments = new ArrayList<>();

    private File blastOut;

    public static void main(String[] args) {
        SwitchBlastParser parser = new SwitchBlastParser(new File("D:\\Users\\Mathilde\\Documents\\KATHRIN" +
                "\\SwitchRegionAnalysis\\PIPELINE_V4\\insertsAnalysis\\blat.out"));
    }

    public SwitchBlastParser(File blastOut){
        this.blastOut= blastOut;
        try {
            parseBlastOutFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void parseBlastOutFile() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(blastOut.getPath()));
        String line = "";
        while ((line = fileReader.readLine()) != null) {
            if (!line.matches("#.*")){
                //System.out.println("Get "+line);
                String[] cells = line.split("\\s+");
                Alignment alignment = new Alignment();
                alignment.setQueryId(cells[0]);
                alignment.setQueryLength(Integer.parseInt(cells[1]));
                alignment.setQueryStart(Integer.parseInt(cells[2]));
                alignment.setQueryEnd(Integer.parseInt(cells[3]));
                alignment.setSubjectId(cells[4]);
                alignment.setSubjectLength(Integer.parseInt(cells[5]));
                alignment.setAlignLength(Integer.parseInt(cells[6]));
                alignment.setPercIdentity(Double.parseDouble(cells[7]));
                alignments.add(alignment);
            }
        }
        fileReader.close();
    }

    public ArrayList<Alignment> getAlignments(){
        return alignments;
    }

    public class Alignment {

        private String queryId;
        private int queryLength;
        private int queryStart;
        private int queryEnd;
        private String subjectId;
        private int subjectLength;
        private int alignLength;
        private double percIdentity;


        public String getQueryId() {
            return queryId;
        }

        public void setQueryId(String queryId) {
            this.queryId = queryId;
        }

        public int getQueryLength() {
            return queryLength;
        }

        public void setQueryLength(int queryLength) {
            this.queryLength = queryLength;
        }

        public int getQueryStart() {
            return queryStart;
        }

        public void setQueryStart(int queryStart) {
            this.queryStart = queryStart;
        }

        public int getQueryEnd() {
            return queryEnd;
        }

        public void setQueryEnd(int queryEnd) {
            this.queryEnd = queryEnd;
        }

        public String getSubjectId() {
            return subjectId;
        }

        public void setSubjectId(String subjectId) {
            this.subjectId = subjectId;
        }

        public int getSubjectLength() {
            return subjectLength;
        }

        public void setSubjectLength(int subjectLength) {
            this.subjectLength = subjectLength;
        }

        public int getAlignLength() {
            return alignLength;
        }

        public void setAlignLength(int alignLength) {
            this.alignLength = alignLength;
        }

        public double getPercIdentity() {
            return percIdentity;
        }

        public void setPercIdentity(double percIdentity) {
            this.percIdentity = percIdentity;
        }

        public Alignment(){

        }
    }
}
