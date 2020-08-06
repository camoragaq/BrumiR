use Data::Dumper;
use Getopt::Std;
use strict;


sub usage {
	print "$0 usage : -a <BCALM.tab> -p <OUT>\n";
	print "Error in use\n";
	exit 1;
}

my %opts = ();
getopts( "a:p:", \%opts );
if ( !defined $opts{a} or !defined $opts{p}) {
	usage;
}



open(FILE,$opts{a}) or die "cannot open file\n";

my $hash=();
my @csvs=();
my $mirs=();
while (my $csv=<FILE>) {
	# body...
	 ($hash,$mirs)=load_exp_by_cmp_miRNA($csv,$hash);
	chomp $csv;
	push(@csvs,$csv);
}

#print Dumper($hash);
my @kmirs=sort keys %{$mirs};
print join(" ",@csvs)."\n";
foreach my $m(@kmirs){
	my @etmp=();
	foreach my $f(@csvs){
		if(defined $hash->{$f}->{$m}){
			my ($bcce)=sort {$hash->{$f}->{$m}->{$b}->{adds} <=> $hash->{$f}->{$m}->{$a}->{adds}} keys %{$hash->{$f}->{$m}};
			#print Dumper($hash->{$f}->{$m});
			#print Dumper($bcce);
			push(@etmp,$hash->{$f}->{$m}->{$bcce}->{adds});
		}else{
			push(@etmp,0);
		}
	}	
	print join(" ",$m,@etmp)."\n";
}







sub load_exp_by_cmp_miRNA{
	my ($file,$hash)=@_;
	chomp $file;
	open(CSV,$file) or die "cannot open CSV file\n";

	while (my $line=<CSV>) {
		# body...
		next if($line=~m/^#/);
		chomp $line;
		my @data=split /\t/,$line;
		####CCID   SCC     UNITIG_ID       MIRBASE KCOV    KC      LENGTH  UNISEQ
#5       4       45964   ssc-miR-181a    1131    24895   22      ACTCACCGACAGCGTTGAATGT
#5       4       36280   ssc-miR-181a    1599    35196   22      CTCACCGACAGCGTTGAATGTT
#5       4       45129   ssc-miR-181a    126     2788    22      CATTCAACGCTGTCGGTGAGTT
		#print Dumper(@data);
		$hash->{$file}->{$data[3]}->{$data[0]}->{adds}+=$data[4];
		if($hash->{$file}->{$data[3]}->{$data[0]}->{maxx}  <$data[4] ){
			$hash->{$file}->{$data[3]}->{$data[0]}->{maxx}=$data[4];
		}
		$mirs->{$data[3]}++;
	}
	#print Dumper($hash);
	return ($hash, $mirs);
}


