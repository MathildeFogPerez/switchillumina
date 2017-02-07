import pysam
import argparse
import sys

parser = argparse.ArgumentParser(
    description='Filter a paired-end BAM file to filter the spanning reads.')
parser.add_argument('input_bamFile', help="input bam file name")
parser.add_argument('insert', help="insert name")

args = parser.parse_args()
if (args.input_bamFile == None or args.insert == None):
    parser.print_help()
    sys.exit()

# Now we use only the 2 extremities of the Switch region
startSwitch = 106050000
endSwitch = 106337000


def isSwitchRegion(chr, start):
    if chr == "chr14" and start >= startSwitch and start <= endSwitch:
        return True;
    return False;


# We read the BAM file where the reads
samfile = pysam.AlignmentFile(args.input_bamFile, "rb")
readsNumber = 0
selectedReads = []
# for each read inside the coordinates we select the spanning reads or the ones inside the insert
# without SA tag
for read in samfile.fetch():
    readsNumber += 1
    # we dont want read with tag XA:Z 2nd alignments
    if not read.has_tag("XA"):
        readId = str(read).strip().split()
        # it has a SA tag in the switch region
        if read.has_tag("SA"):
            saTag = read.get_tag("SA").split(",")
            if isSwitchRegion(saTag[0], int(saTag[1])):
                selectedReads.append(readId[0])
        else:
            selectedReads.append(readId[0])

print "Number TOTAL of reads: ", readsNumber
print "Number of selected inside and spanning reads: ", len(selectedReads)

samfile.close()

# print out the selected regions ids.${COLS[1]}.spanning.txt
f = open("ids." + args.insert + ".inside.spanning.txt", 'w+')
for read in selectedReads:
    print >> f, read

