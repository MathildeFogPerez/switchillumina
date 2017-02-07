# README #


1. Trim the reads in rawData dir
>mkdir CleaningWithTrimGalore
>trim_galore --illumina --paired -q 20 --retain_unpaired --length 99 --length_1 100 --length_2 100 --output_dir ../CleaningWithTrimGalore/ IC-2295_C_208072_G_lib156878_5068_2_1.fastq IC-2295_C_208072_G_lib156878_5068_2_2.fastq &

2.Align the reads to the human genome hg19
Important we DO NOT USE the -M option in order to remove the PCR duplicate reads using Picard tools. This doesnt work with SA flag-> mess everything up
>bwa mem -t 10 /home/13010563@unisi.ch/bwa.kit/hg19.fasta /home/13010563@unisi.ch/KATHRIN/Illumina/C_208072_G/CleaningWithTrimGalore/*_val_1.fq /home/13010563@unisi.ch/KATHRIN/Illumina/C_208072_G/CleaningWithTrimGalore/*_val_2.fq > C_208072_G_notSorted.sam &

>samtools view -Sb  C_208072_G_notSorted.sam > C_208072_G_notSorted.bam
>rm C_208072_G_notSorted.sam
>samtools sort C_208072_G_notSorted.bam C_208072_G

>disown -h jobId

>rm C_208072_G_notSorted.bam 

*****Get some stats
>samtools index C_208072_G.bam
>samtools flagstat *.bam

3. Get the "over covered" regions
>bedtools genomecov -dz -ibam C_208072_G.bam > C_208072_G_depth_v4_bedtools.txt

4. Identify potential inserts among those over covered regions
a.WE FILTER the regions with 2 reads/bp and mini length =50bp (java prog)
>java -jar /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/FindOverCoverRegion.jar C_208072_G 2 50
We keep insert with chimeric reads in 3' and 5' with Switch region and at least 2 mates map in Switch region (ValidateInserts.py)
We dont process regions >= 2000 bps (probably non specific PCR product) and we process only read where map quality >= 5 and not have XA tag (multi mapping reads)
>python /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/ValidateInserts.py C_208072_G.bam C_208072_G 2 &
>disown -h PID

b.WE FILTER the regions with 40 reads/bp and mini length =50bp (java prog)
>java -jar /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/FindOverCoverRegion.jar C_208072_G 40 50
We keep insert with chimeric reads in 3' and 5' with Switch region and at least 2 mates map in Switch region (ValidateInserts.py)
We dont process regions >= 2000 bps (probably non specific PCR product) and we process only read where map quality >= 5 and not have XA tag (multi mapping reads)
>python /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/ValidateInserts.py C_208072_G.bam C_208072_G 40 &
>disown -h PID

c. We merge 1 and 2 and discard the duplicate regions when the coordinates are EXACTLY the same
>java -jar /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/Merge2Beds.jar C_208072_G selectedInsert_C_208072_G_2reads.bed selectedInsert_C_208072_G_40reads.bed

d. We sort the merged file (IMPORTANT: we sort with this function 'sort-bed'!!)
>sort-bed selectedInsert_C_208072_G_merged.bed > selectedInsert_C_208072_G_merged_sorted.bed

e. We map the annotation to our inserts
>bedmap --echo --echo-map-id-uniq --delim '\t' selectedInsert_C_208072_G_merged_sorted.bed /home/13010563@unisi.ch/Homo_sapiens/hg19/gencode/gencode.v19.annotation.exon.gene_shortedV2.bed > selectedInsert_C_208072_G_Annotated.bed

f. To sort the final file:
>sort -k1,1V -k2,2n selectedInsert_C_208072_G_Annotated.bed > selectedInsert_C_208072_G_Annotated_sorted.bed


g. We create a table with all the info: insert coordinates + merge info +  gene/exon names (CalculateInsertCoverage.java)
In bwa.kit/C_208072_G folder:
>java -jar /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/CalculateInsertCoverage.jar C_208072_G selectedInsert_C_208072_G_Annotated_sorted.bed

--->154 inserts 

h. We add an unique insert Id to the tsv file for each insert and we create a 'selectedInsert_C_208072_G_CONTIG.txt' file
that will be used to create contig sequences for each insert (using run_insert.sh bash script).
>java -jar /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/TsvAnnotatedToInsertId.jar C_208072_G


i.We run Trinity for each insert
>/home/13010563@unisi.ch/SCRIPTS/switch/Illumina/./runInserts.sh C_208072_G &


j.We check that we have nice contig for each insert, otherwise we remove them
First we launch a blast for each contig against its insert and the switch region
>/home/13010563@unisi.ch/SCRIPTS/switch/Illumina/AfterTrinitySelectInsert.sh C_208072_G &
Then we parse the blast output file
>java -jar /home/13010563@unisi.ch/SCRIPTS/switch/Illumina/KeepInsertIfContigIsSwitchInsertSwitch.jar C_208072_G
----> 78 inserts !!!!

>rm -f *.insert*


This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact