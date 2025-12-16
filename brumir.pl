
###############################################################################
# Author: Carol Moraga Quinteros 
# Laboratory: ERABLE-UCBL
# Copyright (c)
# year: 2020
###############################################################################

use Data::Dumper;
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/lib";
#local brumir libraries
use UGraph;

sub usage {
	print "$0 usage : -a <fastq-reads> -p <Prefix>\n";
	print "\nAdvanced options:\n";
	print "\t-k kmer size [def:14]\n";
	print "\t-T Tip coverage [def:5]\n";
	print "\t-A Min Abundance BCALM [def:5]\n";
	print "\t-R factor to split unipaths by coverage [def:3]\n";
	print "\t-d Min coverage to report candidates [def:20]\n";
	print "\t-b RFAM 16-mer database [def:rfam-mer.database.txt]\n";
	print "Error in use\n";
	exit 1;
}
#bcalm
#bcalm -in reads.fastq.gz -kmer-size 18 -abundance-min 5

my %opts = ();
getopts( "a:b:k:d:p:T:R:A:", \%opts );
if ( !defined $opts{a} or !defined $opts{p}) {
	usage;
}

if(!defined $opts{k}){
        $opts{k}=14;
}

if(!defined $opts{T}){
	$opts{T}=5;
}

if(!defined $opts{R}){
	$opts{R}=3;
}

my $min_km=20;
if(defined $opts{d}){
	$min_km=$opts{d};
}

my $bcalm_abundance=5;
if(defined $opts{A}){
   $bcalm_abundance = $opts{A};	
}
my $rfam_mers=$Bin."/db/rfam-mer.database.txt";
if(defined $opts{b}){
	$rfam_mers=$opts{b};
}

my $ed_aligner=$Bin."/bin/ed_aligner";
#we check if ed_aligner is accesible
if(-x $ed_aligner == 1){
}else{
	print "We can not find $ed_aligner\n";
	print "Do you compiled? run [make all] at $Bin\n";
	exit 1;
}
my $bcalm=$Bin."/bin/bcalm";
#we check if bcalm is accesible
if(-x $bcalm == 1){
}else{  
        print "We can not find $bcalm\n";
        print "Do you download it? get the BCALM 2.2.2 from: \nhttps://github.com/GATB/bcalm/releases/download/v2.2.2/bcalm-binaries-v2.2.2-Linux.tar.gz\n or \nhttps://github.com/GATB/bcalm/releases/download/v2.2.2/bcalm-binaries-v2.2.2-Mac.tar.gz\nand copy the BCALM program $Bin, we use version 2.2.2\n";
        exit 1; 
}
#we run bcalm
my $cmd="$bcalm -abundance-min $bcalm_abundance -in $opts{a} -kmer-size $opts{k} -out $opts{p} 2>$opts{p}.bcalm.err > $opts{p}.bcalm.log";
run_cmd($cmd);
# we catch if there is an error runing bcalm
if(!-s "$opts{p}.unitigs.fa"){
  print STDERR "$opts{p}.unitigs.fa file was not create\n";
  exit 1;	
}

#start the brumir algorithm
my $ugraph=new UGraph("$opts{p}.unitigs.fa") or die "cannot create the Graph from the BCALM output\n";
#remove the tips from the Graph i.e sequencing errors, by default we perform 3 iterations
print "####### REMOVE TIP 1 #########\n";
$ugraph->remove_tips(24,$opts{T});
print "####### REMOVE TIP 2 #########\n";
$ugraph->remove_tips(24,$opts{T});
print "####### REMOVE TIP 3 #########\n";
$ugraph->remove_tips(24,$opts{T});
#number of neiborhood to consider and maximum length of the path within a SCC (Strongly Connected Component)
print "####### BREAK EDGES BY COVERAGE #########\n";
$ugraph->del_neigboords_by_coverage2($opts{R});
#another round of polishing of the graph
print "####### REMOVE TIP 4 #########\n";
$ugraph->remove_tips(24,$opts{T});
print "####### REMOVE TIP 5 #########\n";
$ugraph->remove_tips(24,$opts{T});
print "####### REMOVE TIP 6 #########\n";
$ugraph->remove_tips(24,$opts{T});
#delete unipaths by topology
print "####### DELETE BY TOPOLOGY #########\n";
$ugraph->delete_by_topology($opts{p});
#recover and build the new unipaths after applying the previous steps.
print "####### ASSEMBLE UNIPATHS #########\n";
$ugraph->assemble_unipaths($opts{k},$opts{p});
#we compute the overlap among the candidates exact match of K-2 bases
print "####### EXACT OVERLAP #########\n";
$ugraph->compute_overlap_candidates_exact(16,$opts{p});#ideally k-2 bases
#we compute the overlap among the candidates allowing a maximum edit distance of 2
print "####### APROXIMATE OVERLAP #########\n";
$ugraph->compute_overlap_candidates_edit(2,$opts{p},$ed_aligner);
#we compute the overlap with the assembled long sequences (L>=40bp)
print "####### OVERLAP TO ASSEMBLED LONG SEQUENCES #########\n";
$ugraph->overlap_to_othersequences(16,$opts{p});#we map mid-low coverage potential unipaths to long sequences
#we match the RFAM kmers(k16 and Freq >= 5 a total of 6523751 kmers)
print "####### OVERLAP TO RFAM DATABASE #########\n";
$ugraph->overlap_candidates_rfam(16,$rfam_mers,$opts{p});
#we write the output the value is the minimal KM value

print "####### OUTPUT miRNA CANDIDATES #########\n";
$ugraph->write_potencial_miRNA($opts{p},$min_km);
print "####### SAVING FINAL GRAPH #########\n";
#function that print the graph in BCALM format
$ugraph->print_graph_bcalm($opts{p});


#function that run a command
sub run_cmd{
  my ($cmd)=@_;
  print $cmd."\n";
  system($cmd); #or die "cannot run the makefile\n";
  if ($? == -1) {
     print "failed to execute: $!\n";
     exit(1);
 }
 elsif ($? & 127) {
     printf "child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
 exit(1);
 }

}

#############################################################
