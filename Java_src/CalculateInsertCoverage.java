package ch.irb.bamtools;

import ch.irb.utils.Consts;

import java.io.*;
import java.math.BigDecimal;
import java.util.ArrayList;

/**
 * Created by Mathilde on 20.07.2016.
 * This class will take in input a bed file with/without the gene/exon info and calculate the mean bp coverage
 * for each region using the donor_depht_v3_bedtools.txt file generated previously by bedtools genomecov function.
 * It can also take into consideration the '*' that tags "duplicated" regions. NOT TRUE ANYMORE
 */
public class CalculateInsertCoverage {
    private static String donor;//= "donor5";
    private static File bedFile;
    //private File bedFile = new File(foldPath + "selectedInsert_" + donor + "_Annotated_sorted.bed");
    //TODO to check contamination:
    //private File bedFile = new File(foldPath + "bedmap_out_3donors_annotated_sorted.bed");

    private String foldPath = System.getProperty("user.dir");
    ;// = "D:\\Users\\Mathilde\\Documents\\workspace2\\BamTools\\data" + Consts.fs;
    private File depthFile;// = new File(foldPath + donor + "_depth_v3_bedtools.txt");
    private ArrayList<Region> regions = new ArrayList<>();

    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("ERROR: must have 2 arguments for the donor and the bed file (_Annotated_sorted.bed)!");
            System.exit(-1);
        }
        donor = args[0];
        bedFile = new File(args[1]);
        CalculateInsertCoverage cal = new CalculateInsertCoverage();
    }

    public CalculateInsertCoverage() {
        try {
            parseBedFile();
            File outFile = new File(foldPath + Consts.fs + "selectedInsert_" + donor + "_bpCoverage_annotated.tsv");
            BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
            out.write("insert\tlength (bp)\tmean bp coverage\tGene name(s)" + Consts.ls);
            for (Region region : regions) {
                double mean = getCoverage(region);
                out.write(region.getName());
                out.write("\t" + region.getLength() + "\t" + mean + "\t");
                if (region.getGeneId() != null) {
                    out.write(region.getGeneId() + Consts.ls);
                } else {
                    out.write(Consts.ls);
                }
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void parseBedFile() throws IOException {
        depthFile = new File(foldPath + Consts.fs + donor + "_depth_v4_bedtools.txt");
        BufferedReader fileReader = new BufferedReader(new FileReader(bedFile));
        String line = "";
        while ((line = fileReader.readLine()) != null) {
            //System.out.println(":"+line);
            String[] cells = line.split("\\t");
            Region region = new Region(cells[0], Integer.valueOf(cells[1]).intValue(), Integer.valueOf(cells[2]).intValue());
            int i = 3;
            if (cells.length > 3) {//there is a gene id
                if (cells[i].contains(";")) {
                    ArrayList<String> selectedGeneExon = new ArrayList<>();
                    String[] geneExonData = cells[i].split(";");
                    ArrayList<String> geneHasExon = new ArrayList<>();
                    for (String data : geneExonData) {
                        if (data.contains("exon")) {
                            String gene = data.split("\\.exon")[0];
                            selectedGeneExon.add(data.replaceAll("\\.exon", "_exon"));
                            geneHasExon.add(gene);
                        }
                    }
                    for (String data : geneExonData) {
                        if (!data.contains("exon") && !geneHasExon.contains(data)) {
                            selectedGeneExon.add(data);
                        }
                    }
                    String geneId = "";
                    int ind = 1;
                    for (String selected : selectedGeneExon) {
                        if (ind < selectedGeneExon.size()) {
                            geneId += selected + ", ";
                        } else {
                            geneId += selected;
                        }
                        ind++;
                    }
                    region.setGeneId(geneId);
                } else {
                    region.setGeneId(cells[i]);
                }
            }

            regions.add(region);
        }
        fileReader.close();
        System.out.println("Number total of regions " + regions.size());
    }

    private double getCoverage(Region region) throws IOException {
        double length = 0;
        double totReads = 0;
        BufferedReader fileReader = new BufferedReader(new FileReader(depthFile));
        String line = "";
        while ((line = fileReader.readLine()) != null) {
            String[] cells = line.split("\t");
            String chr = cells[0];
            if (chr.equals(region.getChr())) {
                int coord = Integer.valueOf(cells[1]).intValue();
                if (coord >= region.getStart() && coord <= region.getEnd()) {
                    BigDecimal bigDec = new BigDecimal(cells[2]);
                    int reads = bigDec.intValue();
                    length += 1;
                    totReads += reads;
                }
            }
        }
        fileReader.close();
        System.out.println("For region " + region.getName() + ", totReads " + totReads + " over length " + length);
        double mean = totReads / length;
        return mean;
    }


}
