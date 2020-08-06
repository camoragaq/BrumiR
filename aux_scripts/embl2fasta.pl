
###############################################################################
# Author: Carol Moraga
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2017
###############################################################################
use Data::Dumper;
use Getopt::Std;
use Bio::SeqIO;
use strict;

sub usage {
   print "$0 usage : -a <emblfile> \n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}


my $in=Bio::SeqIO->new(-format=>'embl',-file=>$opts{a}) or die "cannot open seq file $opts{a}\n";
while(my $seq=$in->next_seq()){
	print ">".$seq->display_id()."\n".$seq->seq()."\n";
}


