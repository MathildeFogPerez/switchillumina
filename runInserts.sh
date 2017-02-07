#!/bin/bash
DONOR=$1 #argument passed in command line
ILLUMINASCRIPTSDIR="/home/userName/IlluminaScriptsDir"

while IFS= read -r line; do
	COLS=() 
	for val in $line ; do
        	COLS+=("$val")
	done
	echo "--> Processing ${COLS[1]}";	

	#1.Extract reads that are inside insert region, pair (proper or not), AND that are not SA (but can have SA), not mate unmmaped with mapQ>= 5
	#=> inside reads and spanning mate pairs in insert chromosome
	samtools view -h -f1 -F2060 -@8 -q5 $DONOR.bam ${COLS[0]} |samtools view -b -h -> $DONOR.${COLS[1]}.inside.spanning.bam
	samtools index $DONOR.${COLS[1]}.inside.spanning.bam
	#Now we call python program to retrieve the spanning and inside reads
	python $ILLUMINASCRIPTSDIR/GetSpanningReads.py $DONOR.${COLS[1]}.inside.spanning.bam ${COLS[1]}

	#2. spanning mate pairs on switch region chromosome with proper paired reads and are SA reads
	samtools view -h -f2051 -F12 -@8 -q5 $DONOR.bam ${COLS[0]} |samtools view -b -h -> $DONOR.${COLS[1]}.spanning2.bam
	samtools view -@8 -h $DONOR.${COLS[1]}.spanning2.bam | grep 'SA:Z:chr14'> $DONOR.${COLS[1]}.spanning2.SAchr14.sam
	grep -v 'XA:' $DONOR.${COLS[1]}.spanning2.SAchr14.sam > $DONOR.${COLS[1]}.spanning2.SAchr14.unique.sam
	cut -f1 $DONOR.${COLS[1]}.spanning2.SAchr14.unique.sam > ids.${COLS[1]}.spanning2.txt

	#3.HERE we use the python script to get the complete reads (not hard-clipped) from the fastq file with this ids lists
	echo "Retrieving the reads in the fastQ file... number total of IDS:"
	cat ids.${COLS[1]}.inside.spanning.txt ids.${COLS[1]}.spanning2.txt > ids.${COLS[1]}.txt
	wc -l ids.${COLS[1]}.txt
	time cat CleaningWithTrimGalore/*_val_1.fq | python $ILLUMINASCRIPTSDIR/extractReadsFromFastQFile.py ids.${COLS[1]}.txt 1 > $DONOR.${COLS[1]}.1.fastq
	time cat CleaningWithTrimGalore/*_val_2.fq | python $ILLUMINASCRIPTSDIR/extractReadsFromFastQFile.py ids.${COLS[1]}.txt 2 > $DONOR.${COLS[1]}.2.fastq
	
	#4. Trinity assembly
	echo "Running Trinity..."
	Trinity --seqType fq --left $DONOR.${COLS[1]}.1.fastq --right $DONOR.${COLS[1]}.2.fastq --CPU 8 --max_memory 20G --normalize_reads --output $DONOR_${COLS[1]}_trinity &

done <"selectedInsert_"$DONOR"_CONTIG.txt"


