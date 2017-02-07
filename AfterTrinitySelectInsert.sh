#!/bin/bash
DONOR=$1 #argument passed in command line
HG19DIR="/home/userName/hg19Dir"
ILLUMINASCRIPTSDIR="/home/userName/IlluminaScriptsDir"

#AFTER runnin run_inserts.sh to make contig via Trinity, we check out the contig launching BLAST
#against the switch region and the consensus insert sequence

echo "Before getting BAM file"
#To get the consensus insert sequences
samtools view -h -L $DONOR.insertListForBlast.bed $DONOR.bam |samtools view -b -h -> $DONOR.inserts.bam
samtools mpileup -d8000 -uf $HG19DIR/hg19.fasta $DONOR.inserts.bam | bcftools call -c - >file.vcf
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz

#for each insert
while IFS= read -r line; do
	COLS=() 
	for val in $line ; do
        	COLS+=("$val")
	done	
	echo "--> Processing ${COLS[1]} with ${COLS[0]}";
	#check that Trinity.fasta file exists, otherwise go to the next insert
	if [ -f ${COLS[1]}_trinity/Trinity.fasta ]; then
		samtools faidx $HG19DIR/hg19.fasta ${COLS[0]} | bcftools consensus file.vcf.gz -o ${COLS[1]}.fa
		#we move the consensus insert fasta file in the trinity output dir
		mv ${COLS[1]}.fa ${COLS[1]}_trinity/
		#we create the subject file
		cat $ILLUMINASCRIPTSDIR/switchRegion_hg19.fasta ${COLS[1]}_trinity/${COLS[1]}.fa > ${COLS[1]}_trinity/switchRegionsInsert.fasta
		#We launch BLAST insert against switch
		blastn -query ${COLS[1]}_trinity/${COLS[1]}.fa -task megablast -dust no -subject $ILLUMINASCRIPTSDIR/switchRegion_hg19.fasta -outfmt "7 qseqid qlen qstart qend sseqid slen length pident evalue bitscore" -out ${COLS[1]}_trinity/insert_blast.out -perc_identity 80
		#We launch BLAST contig against switch + insert
		blastn -query ${COLS[1]}_trinity/Trinity.fasta -task blastn -dust no -subject ${COLS[1]}_trinity/switchRegionsInsert.fasta -outfmt "7 qseqid qlen qstart qend sseqid slen length pident evalue bitscore" -out ${COLS[1]}_trinity/blast.out -perc_identity 70
	fi
done < selectedInsert_"$DONOR"_CONTIG.txt

#we delete the rest of the files
rm -f file.vcf.gz file.vcf.gz.tbi file.vcf $DONOR.inserts.bam

