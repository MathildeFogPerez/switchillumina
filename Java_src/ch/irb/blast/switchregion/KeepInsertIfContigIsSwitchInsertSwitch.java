package ch.irb.blast.switchregion;

import ch.irb.ManageFastaFiles.FastaFileMaker;
import ch.irb.ManageFastaFiles.FastaFileParser;
import ch.irb.utils.Consts;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 Copyright 2017 - Mathilde Foglierini Perez
 This code is distributed open source under the terms of the GNU Free Documention License.

 * This class will give in output the final list of insert for each sample.
 * In input it needs:
 * - the tsv file containing all the related info for each insert : 'selectedInsert_donor_bpCoverage_annotated_forAmiGo.tsv'
 * - the blast output file of the insert against the switch region
 * -the blast output of all contig sequence against the switch region + insert sequences
 * For each insert it will parse those  BLAST files and keep only insert if:
 * - insert is not homologous to the switch
 * - contig contains complete sequence of insert
 * - contig has 50 bp before and after insert that map to switch region
 * The shortest contig that fulfill those criteria will be kept for further analysis.
 * In output it will create a _FINAL.tsv file with the selected inserts + a fasta file with the sequences of the contig
 * and the insert for each insert.
 */
public class KeepInsertIfContigIsSwitchInsertSwitch {

    private static String donor;
    private File tsvFile;
    private ArrayList<Insert> inserts = new ArrayList<>();
    private HashMap<String, String> insertToLine = new HashMap<>();
    private ArrayList<Insert> selectedInserts = new ArrayList<>();
    private ArrayList<String> lines = new ArrayList<>();
    private LinkedHashMap<String, String> selectedIdToSeq = new LinkedHashMap<>();
    private LinkedHashMap<String, String> selectedContigIdToSeq = new LinkedHashMap<>();

    public static void main(String[] args) {
        if (args.length != 1) {
            System.out.println("ERROR: must have 1 arguments for the donor");
            System.exit(-1);
        }
        donor = args[0];
        KeepInsertIfContigIsSwitchInsertSwitch keepInsertAfterContigBlast = new KeepInsertIfContigIsSwitchInsertSwitch();
    }

    public KeepInsertIfContigIsSwitchInsertSwitch() {
        tsvFile = new File("selectedInsert_" + donor + "_bpCoverage_annotated_forAmiGo.tsv");
        try {
            parseTsvFile();
            checkInserts();
            writeOutputTsvFile();
            writeFastaFile();
            writeContigFastaFile();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /*
    The tsv file contains all the related info for each insert
     */
    private void parseTsvFile() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(tsvFile.getPath()));
        String line = "";
        int lineIndex = 0;
        while ((line = fileReader.readLine()) != null) {
            if (lineIndex > 0) {
                String[] cells = line.split("\t");
                String insertId = cells[5];
                lines.add(line);
                insertToLine.put(insertId, line);
                String coord = cells[0];
                Insert insert = new Insert(insertId, coord);
                inserts.add(insert);
                //System.out.println("Processing "+insertId);
            } else {
                lines.add(line);
            }
            lineIndex++;
        }
        fileReader.close();
    }


    private void checkInserts() throws IOException {
        for (Insert insert : inserts) {
            //first we check that the BLAST insert against the Switch doesnt have a HIT
            if (isInsertSequenceNotSwitch(insert)) {
                if (keepThisInsert(insert)) {
                    selectedInserts.add(insert);
                }
            }
            else{
                System.out.println("After BLAST parsing, we remove " + insert.getInsertId());
            }
        }
        System.out.println("Selected inserts: " + selectedInserts.size());
    }

    /*
    We parse the BLAST output of insert against Switch to discard insert homologous to the Switch
     */
    private boolean isInsertSequenceNotSwitch(Insert insert) throws IOException {
        HashMap<String, Query> queryIdToQuery = new HashMap<>();
        File trinityDir = new File(insert.getInsertId() + "_trinity");
        File blastOut = new File(trinityDir.getPath() + Consts.fs + "insert_blast.out");
        //first check we have a blast output for this insert, not the case if there was no contig created by Trinity
        if (!blastOut.exists()) {
            return false;
        }
        //check the blast output for this insert
        SwitchBlastParser switchBlastParser = new SwitchBlastParser(blastOut);
        ArrayList<SwitchBlastParser.Alignment> alignments = switchBlastParser.getAlignments();
        for (SwitchBlastParser.Alignment align : alignments) {
            int miniAlignment = align.getQueryLength()/2;
            //there is a hit with Switch, we dont keep this insert
            if (align.getAlignLength()>= miniAlignment && align.getPercIdentity()>= 80){
                return false;
            }
        }
        return true;
    }

    /*
    This method will parse the BLAST output of all Trinity contigs (for a given insert) against the Switch region and
    the insert consensus sequence. It returns true if we have Switch(50bp)+insert+Switch(50bp)
     */
    private boolean keepThisInsert(Insert insert) throws IOException {
        HashMap<String,Query> queryIdToQuery = new HashMap<>();
        File trinityDir = new File(insert.getInsertId() + "_trinity");
        File blastOut = new File(trinityDir.getPath() + Consts.fs + "blast.out");
        //first check we have a blast output for this insert, not the case if there was no contig created by Trinity
        if (!blastOut.exists()){
            return false;
        }
        //check the blast output for this contig
        SwitchBlastParser switchBlastParser = new SwitchBlastParser(blastOut);
        ArrayList<SwitchBlastParser.Alignment> alignments = switchBlastParser.getAlignments();

        for (SwitchBlastParser.Alignment align : alignments) {
            //we check that the insert is complete inside the contig +/- 1%
            if (align.getSubjectId().contains("chr") && align.getPercIdentity() >= 98) { //the contig match with the insert at 98% identity
                int marging = (int) (align.getSubjectLength() * 0.01);
                int min = align.getSubjectLength() - marging;
                if (align.getAlignLength() <= align.getSubjectLength() && align.getAlignLength() >= min) {
                    int startLength = align.getQueryStart();
                    int endLength = align.getQueryLength() - align.getQueryEnd();
                    if (align.getQueryEnd() < align.getQueryStart()) { //the insert is strand -
                        startLength = align.getQueryLength() - align.getQueryStart();
                        endLength = align.getQueryEnd();
                    }
                    //the flanking regions of the switch are >= 50bp
                    if (startLength >= 50 && endLength >= 50) {
                        String queryId = align.getQueryId() + " len=" + align.getQueryLength();
                        Query query= new Query(queryId);
                        queryIdToQuery.put(queryId,query);
                        int startInsert = align.getQueryStart();
                        int endInsert = align.getQueryEnd();
                        if (align.getQueryEnd() < align.getQueryStart()) { //the insert is strand -
                            startInsert = align.getQueryEnd();
                            endInsert = align.getQueryStart();
                            query.setInsertIsNegativeStrand(true);
                        }
                        query.setStartInsert(startInsert);
                        query.setEndInsert(endInsert);
                    }
                }
            }
        }

        //then we check that we have a match with the switch in the flanking regions
        if (queryIdToQuery.size()>0) {//at least 1 query has an insert in the middle of the contig
            for (SwitchBlastParser.Alignment align : alignments) {
                String id = align.getQueryId() + " len=" + align.getQueryLength();
                //we check that the insert is complete inside the contig +/- 1%
                if (!align.getSubjectId().contains("chr") && queryIdToQuery.containsKey(id)) { //the contig match with the switch
                    Query query = queryIdToQuery.get(id);
                    int start = align.getQueryStart();
                    int end = align.getQueryEnd();
                    if (align.getQueryEnd() < align.getQueryStart()) {
                        start = align.getQueryEnd();
                        end = align.getQueryStart();
                    }
                    int startInsert = query.getStartInsert();
                    int endInsert = query.getEndInsert();
                    int internalMarge=5;
                    if (end <= startInsert + internalMarge) { //left part of the insert
                        if (align.getAlignLength() >= 50) {
                            query.setLeftPartMatchesWithSwitch(true);
                        }
                    } else if (start + internalMarge >= endInsert) { //right part of the insert
                        if (align.getAlignLength() >= 50) {
                            query.setRightPartMatchesWithSwitch(true);
                        }
                    }
                }
            }
        }

        int min=10000;
        //Keep the shortest contig sequence if we have several which fulfill criteria
        Query selectedQuery=null;
        for(String queryId: queryIdToQuery.keySet()){
            Query query = queryIdToQuery.get(queryId);
            if (query.isLeftPartMatchesWithSwitch() && query.isRightPartMatchesWithSwitch() && query.getQueryLength()<min){
                selectedQuery = query;
            }
        }

        if (selectedQuery != null){//we store the information for the insert
            File insertFile = new File(insert.getInsertId() + "_trinity" + Consts.fs + insert.getInsertId() + ".fa");
            FastaFileParser parser = new FastaFileParser(insertFile);
            HashMap<String, String> idToSeq = parser.getFastaIdToSequence();
            String insertSeq = idToSeq.get(insert.getCoord());
            //System.out.println("GET insert seq "+insertSeq+" for coord "+insert.getCoord());
            selectedIdToSeq.put(insert.getInsertId(), insertSeq);
            insert.setInsertSequence(insertSeq);
            insert.setNegativeStrand(selectedQuery.isInsertIsNegativeStrand());
            //Then the contig sequence
            File contigFile = new File(insert.getInsertId() + "_trinity" + Consts.fs + "Trinity.fasta");
            parser = new FastaFileParser(contigFile);
            idToSeq = parser.getFastaIdToSequence();
            String contigSeq = null;
            for (String id : idToSeq.keySet()) {
                String seq = idToSeq.get(id);
                if (id.contains(selectedQuery.getQueryId())) {
                    contigSeq = seq;
                    break;
                }
            }
            selectedContigIdToSeq.put("contig_" + insert.getInsertId(), contigSeq);
            selectedIdToSeq.put("contig_" + insert.getInsertId(), contigSeq);
            insert.setContigSequence(contigSeq);
            System.out.println("For " + insert.getInsertId() + " we keep " + selectedQuery.getQueryId());
            return true;
        }
        return false;
    }

    private void writeOutputTsvFile() throws IOException {
        File outFile = new File(tsvFile.getName().replace(".tsv", "_FINAL.tsv"));
        BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
        out.write(lines.get(0) + Consts.ls);
        for (Insert insert : selectedInserts) {
            String line = insertToLine.get(insert.getInsertId());
            out.write(line + Consts.ls);
        }
        out.close();
    }

    private void writeFastaFile() throws IOException {
        File fastaFile = new File(donor + "_inserts_contigs.fasta");
        FastaFileMaker maker = new FastaFileMaker(fastaFile.getPath(), selectedIdToSeq);
    }

    private void writeContigFastaFile() throws IOException {
        File fastaFile = new File(donor + "_contigs.fasta");
        FastaFileMaker maker = new FastaFileMaker(fastaFile.getPath(), selectedContigIdToSeq);
    }

    private class Query {
        private String queryId;
        private int queryLength;
        private int startInsert;
        private int endInsert;
        private boolean leftPartMatchesWithSwitch = false;
        private boolean rightPartMatchesWithSwitch = false;
        private boolean insertIsNegativeStrand = false;

        public Query (String queryId){
            setQueryId(queryId);
            setQueryLength(Integer.parseInt(queryId.split("len=")[1]));
        }
        public String getQueryId() {
            return queryId;
        }

        public void setQueryId(String queryId) {
            this.queryId = queryId;
        }
        public boolean isInsertIsNegativeStrand() {
            return insertIsNegativeStrand;
        }

        public void setInsertIsNegativeStrand(boolean insertIsNegativeStrand) {
            this.insertIsNegativeStrand = insertIsNegativeStrand;
        }

        public int getQueryLength() {
            return queryLength;
        }

        public void setQueryLength(int queryLength) {
            this.queryLength = queryLength;
        }

        public int getStartInsert() {
            return startInsert;
        }

        public void setStartInsert(int startInsert) {
            this.startInsert = startInsert;
        }

        public int getEndInsert() {
            return endInsert;
        }

        public void setEndInsert(int endInsert) {
            this.endInsert = endInsert;
        }

        public boolean isLeftPartMatchesWithSwitch() {
            return leftPartMatchesWithSwitch;
        }

        public void setLeftPartMatchesWithSwitch(boolean leftPartMatchesWithSwitch) {
            this.leftPartMatchesWithSwitch = leftPartMatchesWithSwitch;
        }

        public boolean isRightPartMatchesWithSwitch() {
            return rightPartMatchesWithSwitch;
        }

        public void setRightPartMatchesWithSwitch(boolean rightPartMatchesWithSwitch) {
            this.rightPartMatchesWithSwitch = rightPartMatchesWithSwitch;
        }


    }


    private class Insert {


        private String insertId;
        private String coord;
        private String insertSequence;
        private String contigSequence;
        private boolean isNegativeStrand = false;

        public Insert(String insertId, String coord) {
            this.insertId = insertId;
            this.coord = coord;
        }

        public String getInsertSequence() {
            return insertSequence;
        }

        public void setInsertSequence(String insertSequence) {
            this.insertSequence = insertSequence;
        }

        public String getContigSequence() {
            return contigSequence;
        }

        public void setContigSequence(String contigSequence) {
            this.contigSequence = contigSequence;
        }

        public boolean isNegativeStrand() {
            return isNegativeStrand;
        }

        public void setNegativeStrand(boolean negativeStrand) {
            isNegativeStrand = negativeStrand;
        }

        public String getInsertId() {
            return insertId;
        }

        public String getCoord() {
            return coord;
        }
    }
}
