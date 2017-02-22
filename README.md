# Switch region PIPELINE (Illumina) #

Copyright (C) 2017  Mathilde Foglierini Perez / Kathrin Pieper

email: mathilde.perez@irb.usi.ch kathrin.pieper@irb.usi.ch

### SUMMARY ###

We have made available here a series of scripts to analyze the Switch region (IGH locus) to find potential DNA insertions.

The scripts are primarily intended as reference for manuscript REF_TO_PUT_HERE rather than a stand-alone application.

The input of the pipeline is 300 bp paired-end reads coming from a target amplicon of the switch region. Data can be found at SRA_REF_TO_PUT_HERE accession number.

These scripts were run on Linux machines. The general overview of the pipeline is represented in 'SuppFig_Pipeline_Illumina.pdf' file.


### LICENSES ###

This code is distributed open source under the terms of the GNU Free Documention License.


### INSTALL ###

Before the pipeline can be run, the following software are required:

a) Python 2.7 https://www.python.org/download/releases/2.7/

b) Pysam https://github.com/pysam-developers/pysam

c) Java JDK 8 https://docs.oracle.com/javase/8/docs/technotes/guides/install/install_overview.html#A1097144

d) FastQC and Trim Galore! http://www.bioinformatics.babraham.ac.uk/projects/index.html

e) Burrows-Wheeler Aligner http://bio-bwa.sourceforge.net/

f) samtools v1.3.1 http://samtools.sourceforge.net/

g) Bedtools v2.26 http://bedtools.readthedocs.io/en/latest/index.html#

h) BEDOPS v2.4.20 https://bedops.readthedocs.io/en/latest/

i) Trinity v2.3.2 https://github.com/trinityrnaseq/trinityrnaseq/wiki

j) BLAST+ v2.5.0 ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/

### PIPELINE ###

First we create a directory for each donor ($DONOR) and move in the raw fastq files.


1. Trim the reads

        $ mkdir CleaningWithTrimGalore

        $ trim_galore --illumina --paired -q 20 --length 99 --output_dir CleaningWithTrimGalore/ raw_1.fastq raw_2.fastq


2. Align the reads to the human genome hg19 (allowing soft-clipped reads)

        $ bwa mem -t 10 /pathToHg19/hg19.fasta CleaningWithTrimGalore/*_val_1.fq CleaningWithTrimGalore/*_val_2.fq | samtools view -bSu - > $DONOR.notSorted.bam

        $ samtools sort $DONOR.notSorted.bam $DONOR

    Get some statistics

        $ samtools index $DONOR.bam
        $ samtools flagstat $DONOR.bam


3. Get the "over covered" regions

        $ bedtools genomecov -dz -ibam $DONOR.bam > $DONOR_depth_v4_bedtools.txt


4. Identify potential inserts among those over covered regions

    *  We filter in the regions mapped by minimum 2 reads/bp and mini length =50bp

             $ java -jar /pathToSwitchIlluminaScripts/FindOverCoverRegion.jar $DONOR 2 50

         We keep insert with chimeric reads in 3' and 5' with Switch region and at least 2 mates map in Switch region

         We discard regions >= 2000 bps (probably non specific PCR product) and we process only read where map quality >= 5  and not have XA tag (multi mapping reads)

             $ python /pathToSwitchIlluminaScripts/ValidateInserts.py $DONOR.bam $DONOR 2


    * We filter in the regions mapped by minimum 40 reads/bp and mini length =50bp

             $ java -jar /pathToSwitchIlluminaScripts/FindOverCoverRegion.jar $DONOR 40 50

         We keep insert with chimeric reads in 3' and 5' with Switch region and at least 2 mates map in Switch region

         We discard regions >= 2000 bps (probably non specific PCR product) and we process only read where map quality >= 5  and not have XA tag (multi mapping reads)

             $ python /pathToSwitchIlluminaScripts/ValidateInserts.py $DONOR.bam $DONOR 40



5. We merge the insert coordinates found with minimum 2 reads/bp coverage and minimum 40 reads/bp coverage.

      If two insert coordinates overlap and the distance between the two is equal or below 10 bp we keep the shortest insert,
      else we keep the longest one.

        $ java -jar /pathToSwitchIlluminaScripts/Merge2Beds.jar $DONOR selectedInsert_$DONOR_2reads.bed selectedInsert_$DONOR_40reads.bed


6. We sort the merged file

        $ sort-bed selectedInsert_$DONOR_merged.bed > selectedInsert_$DONOR_merged_sorted.bed

7. We map the annotation to our inserts

        $ bedmap --echo --echo-map-id-uniq --delim '\t' selectedInsert_$DONOR_merged_sorted.bed /pathToSwitchIlluminaScripts/gencode.v19.annotation.exon.gene_shortedV2.bed > selectedInsert_$DONOR_Annotated.bed

8. We sort the final file:

        $ sort -k1,1V -k2,2n selectedInsert_$DONOR_Annotated.bed > selectedInsert_$DONOR_Annotated_sorted.bed

    and we create a tsv file with all the info (insert coordinates, annotation, reads coverage)

        $ java -jar /pathToSwitchIlluminaScripts/CalculateInsertCoverage.jar $DONOR selectedInsert_$DONOR_Annotated_sorted.bed


9. We add an unique insert Id to the tsv file for each insert and we create a 'selectedInsert_$DONOR_CONTIG.txt' file that will be used to create contig sequences for each insert (used in runInsert.sh bash script).

        $ java -jar /pathToSwitchIlluminaScripts/TsvAnnotatedToInsertId.jar $DONOR


10. For each insert, we will extract the encompassing paired reads and spanning paired reads in order

    to create a contig sequence using Trinity software.

        $ /pathToSwitchIlluminaScripts/./runInserts.sh $DONOR


11. We finally check that we have nice contig for each insert, otherwise we remove them

    First we launch a blast for each insert against the switch region to discard the insert sequences that are homologous with the switch region, and then a blast for each contig against its insert and the switch region to keep only the contig that are Switch/Insert/Switch.

        $ /pathToSwitchIlluminaScripts/AfterTrinitySelectInsert.sh $DONOR

    Then we parse the blast output files

        $ java -jar /pathToSwitchIlluminaScripts/KeepInsertIfContigIsSwitchInsertSwitch.jar $DONOR


    Three files will be produced:

    * 'selectedInsert_$DONOR_bpCoverage_annotated_forAmigo_FINAL.tsv' (can be used to add GO terms to the insert list)
    * '$DONOR_contigs.fasta' which contains a contig for each validated insert
    * '$DONOR_inserts_contigs.fasta' which contains the contig and the insert sequences for each validated insert