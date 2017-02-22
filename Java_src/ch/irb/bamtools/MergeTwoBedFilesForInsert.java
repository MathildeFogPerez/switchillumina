package ch.irb.bamtools;

import ch.irb.utils.Consts;

import java.io.*;
import java.util.ArrayList;

/**
 Copyright 2017 - Mathilde Foglierini Perez
 This code is distributed open source under the terms of the GNU Free Documention License.

 This class is used when we want to merge the 2 bed files coming from the 2 workflows: minimum reads coverage = 2 reads
 and minimum reads coverage = 40 reads. We will remove the duplicated inserts and if 2 insert overlap and the difference
 of overlapping is <= 10bp we keep the shortest one, otherwise we keep the longest one.
 */
public class MergeTwoBedFilesForInsert {
    private static String donor;
    private ArrayList<Region> regions1 = new ArrayList<>();
    private ArrayList<Region> regions2 = new ArrayList<>();
    private ArrayList<Region> toDiscard = new ArrayList<>();
    private static File bedFile1;
    private static File bedFile2;

    public static void main(String[] args) {
        if (args.length != 3) {
            System.out.println("Missing 2 files and the donor in arguments!!");
            System.exit(0);
        }
        donor = args[0];
        bedFile1 = new File(args[1]);
        bedFile2 = new File(args[2]);
        MergeTwoBedFilesForInsert merge = new MergeTwoBedFilesForInsert();
    }

    public MergeTwoBedFilesForInsert() {

        System.out.println("Processing donor: " + donor);
        try {
            regions1 = parseBedFile(bedFile1);
            regions2 = parseBedFile(bedFile2);
            processRegions();
            writeMergedFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private ArrayList<Region> parseBedFile(File bedFile) throws IOException {
        ArrayList<Region> regions = new ArrayList<>();
        BufferedReader fileReader = new BufferedReader(new FileReader(bedFile));
        String line = "";
        while ((line = fileReader.readLine()) != null) {
            //System.out.println(":"+line);
            String[] cells = line.split("\\t");
            Region region = new Region(cells[0], Integer.valueOf(cells[1]).intValue(), Integer.valueOf(cells[2]).intValue());
            regions.add(region);
        }
        fileReader.close();
        return regions;
    }

    private void processRegions() {
        for (Region reg1 : regions1) {
            for (Region reg2 : regions2) {
                if (reg1.getName().equals(reg2.getName())) {
                    toDiscard.add(reg2);
                    System.out.println("This region is duplicated " + reg2.getName());
                } else if (reg1.getChr().equals(reg2.getChr())) {
                    //if there are duplicate region, we remove one of them (shortest or longest one, see method)
                    if (areDuplicateRegions(reg1,reg2)){
                        toDiscard.add(removeDuplicateRegions(reg1,reg2));
                    }
                }
            }
        }

        for (Region reg2 : regions2) {
            if (!toDiscard.contains(reg2)) {
                for (Region reg1 : regions1) {
                    if (!toDiscard.contains(reg1)){
                        if (reg1.getChr().equals(reg2.getChr())) {
                            //if there are duplicate region, we remove one of them (shortest or longest one, see method)
                            if (areDuplicateRegions(reg2,reg1)){
                                toDiscard.add(removeDuplicateRegions(reg2,reg1));
                            }
                        }
                    }
                }
            }
        }
    }



    private void writeMergedFile() throws IOException {
        File outFile = new File("selectedInsert_" + donor + "_merged.bed");
        BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
        for (Region region : regions1) {
            if (!toDiscard.contains(region)) {
                out.write(region.getChr() + "\t" + region.getStart() + "\t" + region.getEnd());
                out.write(Consts.ls);
            }
        }
        for (Region region : regions2) {
            if (!toDiscard.contains(region)) {
                out.write(region.getChr() + "\t" + region.getStart() + "\t" + region.getEnd());
                out.write(Consts.ls);
            }
        }
        out.close();
    }

    private boolean areDuplicateRegions(Region reg1, Region reg2){
        if (reg1.getStart() >= reg2.getStart() && reg1.getStart() <= reg2.getEnd()) {
            return true;
        } else if (reg1.getEnd() >= reg2.getStart() && reg1.getEnd() <= reg2.getEnd()) {
            return true;
        }
        return false;
    }

    private Region removeDuplicateRegions (Region reg1, Region reg2){
        int length1 = reg1.getLength();
        int length2= reg2.getLength();
        int diff = Math.abs(length1-length2);
        if (diff<=10){ //10 bp of diff between the 2 regions, we keep the shortest one and remove the longest one
            if (length1<length2){
                reg1.setDuplicate(true);
                return reg2;
            }
            else{
                reg2.setDuplicate(true);
                return reg1;
            }
        }
        else{ // more than 10bp of diff, we keep the longest one and remove the shortest one
            if (length1>length2){
                reg1.setDuplicate(true);
                return reg2;
            }
            else{
                reg2.setDuplicate(true);
                return reg1;
            }
        }
    }
}
