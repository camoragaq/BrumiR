from collections import defaultdict
from Bio import SeqIO


file=open("plants.mirbase.txt","r")
dplants=defaultdict(set)

for line in file:
 line=line.rstrip()
 idp, *rest=line.split("\t")
 dplants[idp]=rest
 #print(idp,rest)


fasta=open("miRbase.mature.plants.fa", 'w')

#print(dplants)
for seq_record in SeqIO.parse("mature.fa", "fasta"):
 idp, *rest =seq_record.id.split("-")
 if(idp in dplants):
  fasta.write('>{}\n{}\n'.format(seq_record.id, str(seq_record.seq))) 
  print(seq_record.id)
  #print(repr(seq_record.seq))
  #print(len(seq_record))
