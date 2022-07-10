use Data::Dumper;
use SeqComplex;
use Getopt::Std;


sub get_Xmer_com{
	my ($x,$hash)=@_;
	my $newhash=();
foreach my $seq (keys %{$hash}){	
   for(my $i=0; $i < length($seq)-$x; $i++){
		$newhash->{substr($seq, $i, $x)}++;
  }
}

 return ($newhash);
}

sub match_Xmer{
	my ($seq,$x,$hash)=@_;
    my $count=0;
   for(my $i=0; $i < length($seq)-$x; $i++){
	my $xmer=substr($seq, $i, $x);
	if(defined $hash->{$xmer}){
		$count+=$hash->{$xmer};
	}
    }
  
    return $count;
  	
}

sub usage {
	die(qq/
	  Usage example:
    	$0  -i brumir.fasta -d mirna.db.fasta  -m rf.model -s sample.name

	$0 options:

	Mandatory options***:
	 -i brumir fasta results
	 -d miRNA database in fasta format
	 -m random forest model
	 -s sample name for output

	\n/);
	exit 1;
}


my %opts = ();
getopts("i:d:m:s:h:",\%opts);
if($opts{h}){
	usage;
}

if(!defined $opts{i} or !defined $opts{m} or !defined $opts{s} or !defined $opts{d}){
	usage;
}


#we read the known miRNAs
open(FILE, $opts{d}) or die "cannot open miRNA data base $opts{d}\n";

while(my $name=<FILE>){
	my $seq=<FILE>;
	chomp $name;
	chomp $seq;
	$seq=~s/U/T/g;
	for(my $i=0; $i < length($seq)-15; $i++){
		$hash->{substr($seq, $i, 15)}++;
	}
}
close(FILE);
#mer databases from know miRNAs
my $mer7=get_Xmer_com(7,$hash);
my $mer6=get_Xmer_com(6,$hash);
my $mer8=get_Xmer_com(8,$hash);


#we load the BrumiR candidates
open(BR,$opts{i}) or die "cannot open Brumir fasta file $opts{i}\n";
my @brcandidates=();
while(my $name=<BR>){
     my $seq=<BR>;
      chomp $name;
      chomp $seq;
	
	#forward
      my ($n)=split(" ",$name);
	$n=~s/>//;
	my $tmp=();
	$tmp->{id}=$n;
	$tmp->{name}=$name;
	$tmp->{seq}=$seq;
	push(@brcandidates,$tmp);

	for(my $i=0; $i < length($seq)-15; $i++){
		#print substr($seq, $i, 15)."\n";
		push(@{$hashbr->{substr($seq, $i, 15)}},$n);
	}
	#reverse seq
	my $revcomp = reverse $seq;
	$revcomp =~ tr/ATGCatgc/TACGtacg/;
	for(my $i=0; $i < length($revcomp)-15; $i++){
		push(@{$hashbr->{substr($revcomp, $i, 15)}},$n);
	}
}

close(BR);

#15 mer database
my @methods = qw/gc gcs cpg cwf ce cm1 cm2 cm3 cm4 cm5 cm6 ct1 ct2 ct3 ct4 ct5 ct6 cl1 cl2 cl3 cl4 cl5 cl6/;

open(M,">".$opts{s}.".br.matrix.txt") or die "cannot create matrix for classification\n";
print M join("\t","mer","type","in_db","mer6","mer7","mer8",@methods)."\n";


foreach my $s (keys %{$hashbr}){
  my %results=runAllMethods($s);
  my @tmp=();
 foreach my $m (@methods) {
                my $val = shift @{ $results{$m} };
		push(@tmp, $val);
        }
	my $m6=match_Xmer($s,6,$mer6);
	my $m7=match_Xmer($s,7,$mer7);
	my $m8=match_Xmer($s,8,$mer8);
	print M join("\t",$s,"BrumiR",join(",",@{$hashbr->{$s}}),$m6,$m7,$m8,,@tmp)."\n";
}

close(M);

# we do classify the brumir candidades 
my $cmd="Rscript Classify_miRNAs_RF.R -i $opts{s}.br.matrix.txt -m $opts{m} -s $opts{s}";
system($cmd) == 0
    or die "system $cmd failed: $?";
# we filter the fasta file with candidates classified as miRNAs

open(FF,$opts{s}."_brumir_predictions_filtered.txt") or die "cannot open filtered Brumir candidates file $opts{s}_brumir_predictions_filtered.txt\n";
my $fhash=();
while(my $line =<FF>){
	chomp $line;
	my (undef,undef,$id,undef)=split("\t",$line);
	$fhash->{$id}++;
}
close(FF);
#print Dumper($fhash);

#we open the file for output Brumir filter
open(RR,">".$opts{s}."_brumir_predictions_filtered.fasta") or die "cannot create output file\n";

foreach my $s(@brcandidates){
	if(defined $fhash->{$s->{id}}){
	  print RR join("\n",$s->{name},$s->{seq})."\n";
	}
}

close(RR);

