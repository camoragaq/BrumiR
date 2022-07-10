use Data::Dumper;
use SeqComplex;


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


open(FILE, $ARGV[0]) or die;

while(my $name=<FILE>){
	my $seq=<FILE>;
	chomp $name;
	chomp $seq;
	$seq=~s/U/T/g;
	#print $name,$seq."\n";
	
	#$hash->{$seq}++;
	for(my $i=0; $i < length($seq)-15; $i++){
		#print substr($seq, $i, 15)."\n";
		$hash->{substr($seq, $i, 15)}++;
	}
	
}

#mer databases from know miRNAs
my $mer7=get_Xmer_com(7,$hash);
my $mer6=get_Xmer_com(6,$hash);
my $mer8=get_Xmer_com(8,$hash);

#15 mer database
my @methods = qw/gc gcs cpg cwf ce cm1 cm2 cm3 cm4 cm5 cm6 ct1 ct2 ct3 ct4 ct5 ct6 cl1 cl2 cl3 cl4 cl5 cl6/;
print join("\t","mer","type","in_db","mer6","mer7","mer8",@methods)."\n";
foreach my $s (keys %{$hash}){
my %results=runAllMethods($s);
  my @tmp=();
 foreach my $m (@methods) {
                my $val = shift @{ $results{$m} };
		push(@tmp, $val);
        }
	my $m6=match_Xmer($s,6,$mer6);
	my $m7=match_Xmer($s,7,$mer7);
	my $m8=match_Xmer($s,8,$mer8);
	print join("\t",$s,"miRNA",$hash->{$s},$m6,$m7,$m8,,@tmp)."\n";
}

### we generate random 15mer sequences
my @dna = qw(A C G T);
sub random_DNA {
    my $length = shift;
    my $seq = '';
    foreach my $n (1..$length) {
        $seq .= $dna[rand(4)]
    }
    return $seq
}

#we generate random sequences
for(my $i=0; $i< 40000; $i++){
  my $s=random_DNA(15);
   $hrand->{$s}=1;		
}

foreach my $s (keys %{$hrand}){
my %results=runAllMethods($s);
  my @tmp=();
 foreach my $m (@methods) {
                my $val = shift @{ $results{$m} };
		push(@tmp, $val);
        }
 my $times=0; 
 if(defined $hash->{$s}){
	$times=$hash->{$s};
  }
	my $m6=match_Xmer($s,6,$mer6);
	my $m7=match_Xmer($s,7,$mer7);
	my $m8=match_Xmer($s,8,$mer8);
	print join("\t",$s,"rand",$times,$m6,,$m7,$m8,@tmp)."\n";
}
