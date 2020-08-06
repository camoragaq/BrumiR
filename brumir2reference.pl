
###############################################################################
# Author: Carol Moraga
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2019
###############################################################################
use Data::Dumper;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/lib";
use Hairpin;
#use RNA::HairpinFigure;
use SSF;
use strict;


sub usage {
   print "$0 usage : -a <brumir-out>  -b <genome.fa> -t<cores> -p <prefix>  -x <plantmode>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:t:p:x:", \%opts );
if ( !defined $opts{a} or !defined $opts{b} or !defined $opts{p}) {
   usage;
}

my $cpus=1;
if(defined $opts{t}){
	$cpus=$opts{t};
}

my $plant=0;
if(defined $opts{x}){
	$plant=1;
}


my $ed_genome_aligner=$Bin."/bin/ed_genome_aligner";
#we check if ed_aligner is accesible
if(-x $ed_genome_aligner == 1){
}else{
	print "We can not find $ed_genome_aligner\n";
	print "Do you compiled? run [make all] at $Bin\n";
	exit 1;
}

my $rnafold=$Bin."/bin/RNAfold";
#we check if ed_aligner is accesible
if(-x $rnafold == 1){
}else{
	print "We can not find $rnafold\n";
	print "Do you download it? get the RNAfold 2.4.9 from:\nhttps://www.tbi.univie.ac.at/RNA/#download\nand copy the RNAfold program $Bin/bin , we use version 2.4.9\n";
	exit 1;
}


#### we run the ed_genome_aligner program
my $cmd="$ed_genome_aligner $opts{a} $opts{b} > $opts{p}.hits.txt";
#activate plant mode
if($plant == 1){
  $cmd="$ed_genome_aligner -p $opts{a} $opts{b} > $opts{p}.hits.txt";
}

run_cmd($cmd);
#### we prepare the file to run the RNAfold program
$cmd="grep \"^hsp\" $opts{p}.hits.txt | awk \'{print \">\"\$1\"\\n\"\$2}\' > $opts{p}.hits.fa";
run_cmd($cmd);

if(-s "$opts{p}.hits.rnafold"){
  print STDERR "$opts{p}.hits.rnafold file exist; deleting it\n";
  $cmd="rm -f $opts{p}.hits.rnafold";
  run_cmd($cmd);
}
### we run the RNAfold program
$cmd="$rnafold -j$cpus  --noPS -i $opts{p}.hits.fa  --outfile=$opts{p}.hits.rnafold";
run_cmd($cmd);
### we parse the RNAfold output
my $fold= new Hairpin("$opts{p}.hits.rnafold","$opts{p}.hits.txt");
##we check every structure and save those that have a plausible hairpin
my $ssf_filter=();
my $total_deleted=0;
open(FILTER,">".$opts{p}.".passfilter.txt") or die "cannot create oputput file\n";
print FILTER join("\t","#hsp","miRNA","chr","start","stop","MFE","B","E", "H", "I","K","M", "S","X","SEGMENTS","Precursor_Seq")."\n";
#hsp	ath-miR156a-5p	1	9503597	9503706	-10.50	1	2	2	2	0	0	5	0	TGAATCTCAAGATATCCACCGTGCTCTCTCTCTTCTGTCAACCTCTTCGGATCCCCTGGCCCAACCACATGTGCAGCCATTTTCTCTACTCTGTTCATATGATGTTGTA

open(NFILTER,">".$opts{p}.".nonpassfilter.txt") or die "cannot create oputput file\n";
print NFILTER join("\t","#hsp","miRNA","chr","start","stop","MFE","B","E", "H", "I","K","M", "S","X","SEGMENTS","Precursor_Seq")."\n";

while(my $f=$fold->next_rnafold_entry()){
    my $isdeleted=0;
    my $st=new SSF($f);
    #we filter the Secondary structure based on the following filters derived from the analysis of precursor sequences present on mirBase
    # Energy
    if(abs($f->{mfe}) >= 15 and abs($f->{mfe}) <= 80){
      #we compute the secondary structure features
      #$st=new SSF($f);
      # paired "Stem"     S
      # Multiloop         M
      # Internal loop     I
      # Bulge             B
      # Hairpin loop      H
      # pseudoKnot        K
      # dangling End      E
      # eXternal loop     X
      #the hash table contains the secondary structure feature counts
      my ($hss)=$st->ssf_hash();

      if($hss->{"H"} eq 1){
          #feature that should be always 0
          if(($hss->{"K"} == 0) and ($hss->{"M"} == 0) and ($hss->{"X"} == 0)){
          #feature that shold be only 1
              if($hss->{"SEGMENTS"} == 1){
                    #we check some feature that shold be less than the Expected
                    if(($hss->{"B"} <= 5) and ($hss->{"E"} <= 3) and ($hss->{"I"} <= 10)) {
                        #is a good candidate.... not deleted using current strcuture filters

                    }else{
                      #deleted by having more than 3 bulges
                      $ssf_filter->{"B"}++ if($hss->{"B"} > 5);
                      #deleted by having more than 2 dangling ends
                      $ssf_filter->{"E"}++ if($hss->{"E"} > 3);
                      #deleted by having more than 5 Internal loops
                      $ssf_filter->{"I"}++ if($hss->{"I"} > 10);
                      $isdeleted=1;
                    }
              }else{
                  $ssf_filter->{"SEGMENTS"}++ if($hss->{"K"} !=1);
                  $isdeleted=1;
              }

          }else{
                #deleted by having a pseudoKnot
                $ssf_filter->{"K"}++ if($hss->{"K"} > 0);
                #deleted by having multiloop
                $ssf_filter->{"M"}++ if($hss->{"M"} > 0);
                #deleted by having an external loop
                $ssf_filter->{"X"}++ if($hss->{"X"} > 0);
                $isdeleted=1;
          }

      }else{
        #deleted by hairpin
        $ssf_filter->{"H"}++;
        $isdeleted=1;
      }

    }else{
      #deleted by low energy
      $ssf_filter->{"MFE"}++;
      $isdeleted=1;
    }
   $total_deleted+=$isdeleted;
  #structure was not deleted;
   if(!$isdeleted){
      #we have to report the structure is a good one
      print FILTER join("\t",split("::",$f->{h}),$f->{mfe}, @{$st->plain_ssf()},$f->{seq})."\n";
   }else{
     print NFILTER join("\t",split("::",$f->{h}),$f->{mfe}, @{$st->plain_ssf()},$f->{seq})."\n";
   }

}
close(FILTER);
print STDERR "Total deleleted candidate hits $total_deleted\n";
foreach my $t(keys %{$ssf_filter}){
  print STDERR join(" ",$t, $ssf_filter->{$t})."\n";
}



#function that run a command
sub run_cmd{
  my ($cmd)=@_;
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
