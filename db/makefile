.DELETE_ON_ERROR:
#makefile to build the 16mer databas from RFAM and miRbase
#external software needed to build the database
#kmc k-mer counter:  https://github.com/refresh-bio/KMC
KMC=/PATH/KMC/bin/kmc
KMCD=/PATH/KMC/bin/kmc_dump
#bioawk 
## git clone https://github.com/lh3/bioawk.git
BAWK=./bioawk/bioawk

download.rfam.txt:
	awk '{print "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/"$$1}' lst.txt > rfam.urls.txt
	wget -c -i rfam.urls.txt
	mkdir -p rfam
	mv *.fa.gz rfam
	touch $@

# the file RFAM-miRNA-FlAM.txt contains all the 529 miRNA families stored in mirbase
rfam.not_miRNA.fa: download.rfam.txt
	ls rfam/*.fa.gz > rfam.lst
	awk '{print $$1}' RFAM-miRNA-FAM.txt > RFAM.miRNA.lst
	grep -v -f RFAM.miRNA.lst rfam.lst > rfam.filter.lst 
	cat rfam.filter.lst | xargs gzip -dc  > $@	

#we consider only k-mers with freq >= 5
rfam.kmc_suf:rfam.not_miRNA.fa
	${KMC} -fm  -n100 -k16 -ci5 $< rfam .
#we dump the k16mers
rfam.not_miRNA.k16f5.txt:rfam.kmc_suf
	${KMCD} rfam $@


#we remove the mature sequences from miRbase
mature_mirbase.kmc_suf:rfam.not_miRNA.k16f5.txt
	wget -c ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
	${KMC} -fm -n100 -k16 -ci1 mature.fa.gz mature_mirbase .

mature_mirbase.k16f1.txt:mature_mirbase.kmc_suf
	${KMCD} mature_mirbase fwd.tmp
	./bioawk/bioawk  '{print revcomp($$1)"\t"$$2}' fwd.tmp > rev.tmp
	cat fwd.tmp rev.tmp > $@
	-rm -f fwd.tmp rev.tmp 
#we remove the kmers in common 
rfam.not_miRNA.k16f5.not_mirbase.txt:rfam.not_miRNA.k16f5.txt mature_mirbase.k16f1.txt
	awk 'BEGIN{while (getline < "mature_mirbase.k16f1.txt"){a[$$1]=$$2}}{if($$1 in a){}else{print $$0}}' rfam.not_miRNA.k16f5.txt > $@

all: download.rfam.txt rfam.not_miRNA.fa  rfam.not_miRNA.k16f5.txt mature_mirbase.kmc_suf mature_mirbase.k16f1.txt rfam.not_miRNA.k16f5.not_mirbase.txt
