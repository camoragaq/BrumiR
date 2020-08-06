#Check that the secondary structure of the miRNA precursor is  Hairpin like.
package Hairpin;
#MirDeep2 script that aim to do the same
#mirdeep2/src/select_for_randfold.pl

use strict;
#create the object
sub new{
  my ($packagename, $rnafold_out, $aligner_out) = @_;
  #print Dumper($unitigs);
  my $self = { rnafold => $rnafold_out,aligner=>$aligner_out};

  bless ($self, $packagename);
  ## create filehandle
  my $filehandle = undef;

	if (ref $rnafold_out eq 'IO::Handle') {
		$filehandle = $rnafold_out;
	}
	else {
		open ($filehandle, $rnafold_out) or die "Error: Couldn't open $rnafold_out\n";
	}
     $self->{fhRNA} = $filehandle;

  return ($self);
}

#### parse this kind of entry
#>hsp_0_0_78861_78971
#GUAAUUUGCGCAAGUAAUUAACAAAAAAAAAUGUAACAAUAAAUGGUGAUCAGAGGCAAUAACAUGUUGGGGGGGAUGAAGCCUGGUCCGAGGAUACUCUCUAUGAUCAC
#((...((((....))))...)).......................(((((((.((((((((...)))))....((((........))))(((....)))))).))))))) (-18.70)
sub next_rnafold_entry{
	my $self=shift;
	my $fh=$self->{fhRNA};
	my $header=<$fh>;
	#there is no more text in the file
	if(!defined $header){
		#we reach the end of the filehandler, must close
		close($self->{fhRNA});
		$self->{fhRNA}=undef;
		return undef;
	}
	my $seq=<$fh>;
	my $smfe=<$fh>;
	my ($struct,$mfe)=split(" ",$smfe);
	if($smfe=~m/\((\s*[-,\s]\d+\.\d+)\)/){
                    $mfe = $1;
        }
	chomp $header;
	$header=~s/>//;
	chomp $seq;
	chomp $mfe;
	#we replace the U->T
	$seq=~s/U/T/g;
	#print join(" ",$header,$seq,$struct,$mfe)."\n";
	my $foldobj= undef;
	$foldobj->{h}=$header;
	$foldobj->{seq}=$seq;
	$foldobj->{struct}=$struct;
	$foldobj->{mfe}=$mfe;
	#we return the foldojb
	return $foldobj;
}

### functions that build and parse the secondary structure




1; #EOM
