#!/usr/bin/python
import sys

#get fname from parameter
idfile=sys.argv[1]
pairIndex=sys.argv[2]
#load ids
ids = set( x.strip() for x in open(idfile) )

#read stdid 
handle=sys.stdin  

while ids:

  #parse fastq
  idline=handle.readline().split(' ')
  seq   =handle.readline()
  spacer=handle.readline()
  quals =handle.readline()
  #check
  completeId = str(idline[0])	
  id=completeId[1:] # remove first character @, everything except the last item :-1
  #print "process id ",id
  if id in ids:
    #print fastq
    toJoin = [completeId,'/',pairIndex,"\n"]
    firstline ="".join(toJoin)
    #print "found id and transform to",firstline
    sys.stdout.write( '%s%s%s%s' % ( firstline, seq, spacer, quals ) )
    #sys.stdout.write( '%s%s%s%s' % ( idline, seq, spacer, quals ) )
    #update ids
    ids.remove( id )

#time cat PL429_q10_filtered_1.fastq | ./parse PL429_q10_filtered_1.fastq.lids > PL429_q10_filtered_1.fastq.lids.fqcpp
