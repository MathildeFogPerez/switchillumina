package ch.irb.bamtools;

import ch.irb.utils.Consts;

import java.io.*;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import static ch.irb.utils.Consts.fs;
import org.apache.commons.math3.stat.descriptive.rank.Median;
/**
 Copyright 2017 - Mathilde Foglierini Perez
 This code is distributed open source under the terms of the GNU Free Documention License.

 This class will parse the depth file produced by bedtools genomecov and select the genomic range with a
 minimum lenght (minLength variable) and a minimum cover (minCover variable) passed in arguments.
 It creates a tsv and a bed file in output with the selected genomic ranges (called "over covered region").
 */

public class FindOverCoveredRegions {
    private static String donor ="MGB47";
    private static int minCover = 2;
    private static int minLength =20;
    private String foldPath= System.getProperty("user.dir")+Consts.fs;
    private File depthFile = new File(foldPath+donor+"_depth_v4_bedtools.txt");


    public static void main(String[] args){
        if (args.length!=3){
            System.out.println("ERROR: must have 3 arguments for the donor, mini read coverage and mini length (bp)!");
            System.exit(-1);
        }
        donor = args[0];
        minCover = Integer.parseInt(args[1]);
        minLength = Integer.parseInt(args[2]);
        try {
            FindOverCoveredRegions find = new FindOverCoveredRegions();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public FindOverCoveredRegions() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(depthFile));
        String line = "";
        double totReads=0;
        double length =0;
        String chrRef=null;
        int start=-1;
        int end=-1;
        ArrayList<Double> allReadsList = new ArrayList<>();
        ArrayList<String> overCoveredRegion = new ArrayList<>();
        LinkedHashMap<String,Integer> regionToReads = new LinkedHashMap<>();
        HashMap<String,Double> regionToMeanReads = new HashMap<>();
        HashMap<String,Integer> regionToLength = new HashMap<>();
        while ((line = fileReader.readLine()) != null) {
            String[] cells = line.split("\t");
            BigDecimal bigDec = new BigDecimal(cells[2]);
            int reads =bigDec.intValue();
            allReadsList.add(new Double(reads));
            if (reads>=minCover) {
                String chr = cells[0];
                int coord = Integer.valueOf(cells[1]).intValue();
                regionToReads.put(chr+"\t"+coord,new Integer(reads));
                if (start==-1){ //we initialize the variable
                    start= coord;
                    end= coord;
                    chrRef = chr;
                    totReads=reads;
                    length=1;
                }
                else {
                    if (chr.equals(chrRef) && coord == end + 1) {
                        end = coord;
                        totReads +=reads;
                        length +=1;
                    }
                    else{//we finished to have a over covered region
                        //we store the data if size is (>= 20bp)
                        if (length>=minLength) {
                            int startToWrite = start+1; //Be careful here we have to add+1 in order to work properly with IGV
                            int endToWrite = end +1; //Be careful here we have to add+1 in order to work properly with IGV
                            String region = chrRef + ":" + startToWrite + "-" + endToWrite;
                            overCoveredRegion.add(region);
                            double mean = totReads /length;
                            regionToMeanReads.put(region,new Double(mean));
                            regionToLength.put(region,Integer.valueOf((int) length));
                        }
                        //we reinitialize start and end
                        start= coord;
                        end= coord;
                        chrRef = chr;
                        totReads=reads;
                        length=1;
                    }
                }
            }
        }
        fileReader.close();
        //we store the last over covered region
        if (length>=minLength) {
            int startToWrite = start+1;
            int endToWrite = end +1;
            String region = chrRef + ":" + startToWrite + "-" + endToWrite;
            overCoveredRegion.add(region);
            double mean = totReads /length;
            regionToMeanReads.put(region,new Double(mean));
            regionToLength.put(region,Integer.valueOf((int) length));
        }
        System.out.println("Number of over covered regions is "+overCoveredRegion.size());
        double[] allReads=new double[allReadsList.size()];
        int i=0;
        for (Double read: allReadsList){
            allReads[i]=read.doubleValue();
            i++;
        }
       // double median = new Median().evaluate(allReads);
        //double mean = allReads/allBase;
        //System.out.println("Number of total reads "+allReads+" for "+allBase+" bases, and median is "+median);

        File outFile = new File(foldPath+"overCoveredRegion_"+donor+"_"+minCover+"readsCoverage.tsv");
        BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
        out.write("# thresholds minCover= "+minCover+" reads/bp, and minLength= "+minLength+" bp"+Consts.ls);
        out.write("Index\tRegion\tSeq length (bp)\tReads mean per seq"+Consts.ls);
        int index=1;
        for (String key: overCoveredRegion){
            out.write(index+"\t"+key);
            out.write("\t"+regionToLength.get(key)+"\t"+regionToMeanReads.get(key).intValue()+Consts.ls);
            index++;
        }
        out.close();

        //we create a bed file that will be used further in the analysis
        File outFile2 = new File(foldPath+"overCoveredRegion_"+donor+"_"+minCover+"reads_V4.bed");
        BufferedWriter out2 = new BufferedWriter(new FileWriter(outFile2));
        for (String key: overCoveredRegion){
            String region = key.replace(":","\t").replace("-","\t");
            out2.write(region+Consts.ls);
            //out2.write(region+"\n"); //TODO Consts.ls or "\n" for LINUX!!
        }
        out2.close();
    }
}
