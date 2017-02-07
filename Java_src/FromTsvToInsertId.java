package ch.irb.bamtools;

import ch.irb.utils.Consts;

import java.io.*;
import java.util.ArrayList;

/**
 * Created by Mathilde on 04.11.2016.
 * This class takes in input the donor/sample name to read the tsv file (i.e.: selectedInsert_2402_G_bpCoverage_annotated.tsv)
 * produced by the pipeline and create a insert id for each insert + a column with a shorten gene id to use for AmiGO
 * It also create a txt file that contains the insert coordinates + insert id used to make the contig
 * and a bed file used to make the insert consensus sequence.
 */
public class FromTsvToInsertId {

    private File tsvFile;
    private static String donor;
    private ArrayList<String> insertsForContig = new ArrayList<>();
    private ArrayList<String> linesToWrite = new ArrayList<>();

    public static void main(String[] args){
        if (args.length != 1) {
            System.out.println("ERROR: must have 1 arguments for the donor!");
            System.exit(-1);
        }
        donor = args[0];
        FromTsvToInsertId fromTsvToInsertId = new FromTsvToInsertId();
    }

    public FromTsvToInsertId(){
        tsvFile = new File("selectedInsert_"+donor+"_bpCoverage_annotated.tsv");
        try {
            parseTsvFile();
            createInsertFileForContig();
            writeBedFile();
            addColumnToTsvFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void parseTsvFile() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(tsvFile.getPath()));
        String line = "";
        int lineIndex=0;
        while ((line = fileReader.readLine()) != null) {
            if (lineIndex > 0) {
                String insertId = "insert"+lineIndex;
                String[] cells = line.split("\t");
                if (cells.length>3){
                    String geneId = getShortGeneId(cells[3]);
                    insertId+= "_"+geneId;
                    linesToWrite.add(line+"\t"+geneId+"\t"+insertId);
                }
                else {
                    linesToWrite.add(line+"\t\t"+insertId);
                }
                insertsForContig.add(cells[0]+"\t"+insertId);
            }
            lineIndex++;
        }
        fileReader.close();
    }

    private void createInsertFileForContig() throws IOException {
        File outFile = new File("selectedInsert_"+donor+"_CONTIG.txt");
        BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
        for (String line: insertsForContig){
            out.write(line+Consts.ls);
        }
        out.close();
    }

    private void writeBedFile()throws IOException{
        File outFile = new File(donor+".insertListForBlast.bed");
        BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
        for (String insert : insertsForContig) {
            String coo = insert.split("\\s")[0];
            String[] coord = coo.split(":");
            out.write(coord[0]+"\t"+coord[1].split("-")[0]+"\t"+coord[1].split("-")[1]+"\t"
                    + insert.split("\\s")[1] + Consts.ls);
        }
        out.close();
    }

    private void addColumnToTsvFile()throws IOException {
        File outFile = new File(tsvFile.getName().replace(".tsv","_forAmiGo.tsv"));
        BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
        out.write("Insert Coordinates\tLength (bp)\tMean bp coverage\tGene name(s)\tGeneId for AmiGO\tInsert id"+Consts.ls);
        for (String line: linesToWrite){
            out.write(line+Consts.ls);
        }
        out.close();
    }

    private String getShortGeneId(String longGeneId){
        String geneId = longGeneId;
        if (longGeneId.contains(",")) {
            geneId = longGeneId.split(",")[0].trim();
            //System.out.println("gene id with , is now "+geneId);
        }
        if (geneId.contains("_")) {
            geneId = geneId.split("_")[0];
        }
        return geneId;
    }
}
