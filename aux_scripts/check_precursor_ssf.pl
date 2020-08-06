
###############################################################################
# Author: Carol Moraga
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2019
###############################################################################
use Data::Dumper;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin";
use Hairpin;
use RNA::HairpinFigure;
use SSF;
use strict;


sub usage {
   print "$0 usage :  -f <rnafold file>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "f:", \%opts );
if ( !defined $opts{f}) {
   usage;
}


### we parse the RNAfold output
my $fold= new Hairpin($opts{f},"bla.hits.txt");
print join(" ","id","MFE","B","E", "H", "I","K","M", "S","X","SEGMENTS")."\n";
##we check every structure and save those that have a plausible hairpin
while(my $f=$fold->next_rnafold_entry()){
    my $st=new SSF($f);
    #method that describe in a table the secondary structure features
    print join(" ",$f->{h},$f->{mfe}, @{$st->plain_ssf()})."\n";
    #print Dumper($st->ssf_hash());
}
