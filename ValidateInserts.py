import pysam
import argparse
import sys

parser = argparse.ArgumentParser(description='Filter a paired-end BAM file to find over cover regions such as mini cover >= a threshold.')
parser.add_argument('input_bamFile', help="input bam file name")
parser.add_argument('donor', help="donor")
parser.add_argument('coverage',type=int, help="minimum coverage")


args = parser.parse_args()
if(args.input_bamFile == None or args.donor == None or args.coverage == None) :
    parser.print_help()
    sys.exit()

#Now we use only the 2 extremities of the Switch region
startSwitch = 106050000
endSwitch = 106337000


def isSwitchRegion(chr, start):
    if chr == "chr14" and start >= startSwitch and start <= endSwitch:
        return True;
    return False;


# We read the BAM file
samfile = pysam.AlignmentFile(args.input_bamFile, "rb")

insertsConfirmed = 0
selectedInserts = []
# First parse the file with the potential insert coordinates
#with open ("nfkbInsert.bed") as f:
with open("overCoveredRegion_" + args.donor + "_" + str(args.coverage) + "reads_V4.bed")as f:  
    index = 0
    for line in f:
        index += 1
        I = line.strip().split()
        startI = int(I[1])
        endI = int(I[2])
        length = endI - startI
        # for each potential insert
        # print "Processing insert ", index, " on ", I[0], " with START ", startI, " and END ", endI
        # check that the "insert" is not in the switch and it is not longer than 2000bp
        if not isSwitchRegion(I[0], startI) and length < 2000:
            hasCondition1 = False
            hasCondition2 = False
            hasCondition3 = False
            read1HasCondi3 = False
            read2HasCondi3 = False
            for read in samfile.fetch(I[0], startI, endI):

                # we go to the next insert if this one is already confirmed
                if hasCondition1 and hasCondition2 and hasCondition3:
                    break

                #the read mapping quality has to be above 5 and we dont have read with tag XA:Z 2nd alignments
                if read.mapping_quality >= 5 and not read.has_tag("XA"):
                    # we check that at least the 5' end has this SA or the 3' end
                    if not hasCondition1 or not hasCondition2:
                        if read.has_tag("SA"):
                            if not hasCondition1:
                                saTag = read.get_tag("SA").split(",")
                                if isSwitchRegion(saTag[0], int(saTag[1])):
                                    startRead = read.reference_start + 1  # to add
                                    if startRead == startI:
					#print "has condition 1: 5' is chimeric"
                                        hasCondition1 = True
                            elif not hasCondition2:
                                saTag = read.get_tag("SA").split(",")
                                if isSwitchRegion(saTag[0], int(saTag[1])):
                                    endRead = read.reference_start + read.query_length
                                    if endRead == endI:
                                        hasCondition2 = True
					#print "has condition 2: 3' is chimeric"

                    # check that at least two mate2 is aligned on the switch region
                    if not hasCondition3:
                        if not read.mate_is_unmapped:
                            mate = samfile.mate(read)
                            if mate.reference_name == "chr14":
                                startMate = mate.reference_start + 1
                                if (isSwitchRegion(mate.reference_name, startMate)):
                                    if not read1HasCondi3:
                                        read1HasCondi3 = True
					#print "has condi 3: 1rst mate is in the switch"
                                    else:  # this is a second read that has its mate in the Switch
                                        hasCondition3 = True
                                        #print "has condi 3: second mate is on switch for read ",read.query_name

            if hasCondition1 and hasCondition2 and hasCondition3:
                insertsConfirmed += 1
                reg = "".join((I[0], "\t", I[1], "\t", I[2]))
                print "->>>>>>>>> insert ", index, ": ", reg, " confirmed!"
                selectedInserts.append(reg)

print "Number of inserts with conditions 1, 2 and 3: ", insertsConfirmed

samfile.close()

# print out the selected regions
f = open("selectedInsert_" + args.donor + "_" + str(args.coverage) + "reads.bed", 'w+')
for insert in selectedInserts:
    print >> f, insert

