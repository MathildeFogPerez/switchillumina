package ch.irb.ManageFastaFiles;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Scanner;
import java.util.regex.Pattern;

/**
 Copyright 2017 - Mathilde Foglierini Perez
 This code is distributed open source under the terms of the GNU Free Documention License.

 FastaParser class, need fasta file in input and return a HashMap idToSeq
 */
public class FastaFileParser {

	private HashMap<String, String> fastaIdToSequence = new HashMap<>();
	private LinkedHashMap<String, String> sameOrderFastaIdToSequence = new LinkedHashMap<>();
	private File file = null;

	public FastaFileParser(File file) throws IOException {
		this.file = file;
		parseFile();
	}

	private void parseFile() throws IOException {
		ArrayList<String> uniqueIds = new ArrayList<>();
		ArrayList<String> sequences = new ArrayList<>();
		Scanner scanner = new Scanner(new FileReader(file));
		scanner.useDelimiter(">");
		// first use a Scanner to get each fasta entry
		while (scanner.hasNext()) {
			int index = 0;
			String fastaId = null;
			String sequence = "";
			Scanner scan = new Scanner(scanner.next());
			scan.useDelimiter(Pattern.compile("([\n]|(\r\n))+"));
			while (scan.hasNext()) {
				String line = scan.next();
				if (index == 0) { // we have the id
					fastaId = line.trim();
				}
				else {
					String seqWithGaps = line.replaceAll("\\s+", "");
					sequence += seqWithGaps; //TODO we dont remove the '-'!!! before .replaceAll("-", "")
				}
				index += 1;
			}
			scan.close();
			// System.out.println("FastaId: "+fastaId+" sequence "+sequence);
			String seq = sequence.toUpperCase().trim();
			if (uniqueIds.contains(fastaId)){
				try {
					throw new Exception("This ID "+fastaId+" is a duplicate!");
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}else{
				uniqueIds.add(fastaId);
			}
			fastaIdToSequence.put(fastaId, seq);
			sameOrderFastaIdToSequence.put(fastaId, seq);
			if (!sequences.contains(seq)) {
				sequences.add(seq);
			}
		}
		scanner.close();
		//System.out.println("Number of fastaIds is: " + fastaIdToSequence.size());
		//System.out.println("Number of unique sequences is: " + sequences.size());
	}

	
	public HashMap<String, String> getFastaIdToSequence(){
		return fastaIdToSequence;
	}
	
	public LinkedHashMap<String, String> getSameOrderFastaIdToSequence(){
		return sameOrderFastaIdToSequence;
	}
	
	public HashMap<String, ArrayList<String>> getSeqToFastaIds(){
		HashMap<String, ArrayList<String>> seqToIds = new HashMap<>();;
		for (String id: fastaIdToSequence.keySet()){
			String seq = fastaIdToSequence.get(id);
			ArrayList<String> ids = new ArrayList<>();
			if (seqToIds.containsKey(seq)){
				ids = seqToIds.get(seq);
			}
			ids.add(id);
			seqToIds.put(seq, ids);
		}
		return seqToIds;
	}
	
}
