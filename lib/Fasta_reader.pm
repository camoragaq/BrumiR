
# lightweight fasta reader capabilities:
package Fasta_reader;

use strict;

sub new {
    my ($packagename, $fastaFile) = @_;

	## note: fastaFile can be a filename or an IO::Handle


    my $self = { fastaFile => undef,,
				 fileHandle => undef };

    bless ($self, $packagename);

    ## create filehandle
    my $filehandle = undef;

	if (ref $fastaFile eq 'IO::Handle') {
		$filehandle = $fastaFile;
	}
	else {

		open ($filehandle, $fastaFile) or die "Error: Couldn't open $fastaFile\n";
		$self->{fastaFile} = $fastaFile;
	}

	$self->{fileHandle} = $filehandle;

    return ($self);
}



#### next() fetches next Sequence object.
sub next {
    my $self = shift;
    my $orig_record_sep = $/;
    $/="\n>";
    my $filehandle = $self->{fileHandle};
    my $next_text_input = <$filehandle>;

	if (defined($next_text_input) && $next_text_input !~ /\w/) {
		## must have been some whitespace at start of fasta file, before first entry.
		## try again:
		$next_text_input = <$filehandle>;
	}

	my $seqobj = undef;

	if ($next_text_input) {
		$next_text_input =~ s/^>|>$//g; #remove trailing > char.
		$next_text_input =~ tr/\t\n\000-\037\177-\377/\t\n/d; #remove cntrl chars
		my ($header, @seqlines) = split (/\n/, $next_text_input);
		my $sequence = join ("", @seqlines);
		$sequence =~ s/\s//g;

		$seqobj = Sequence->new($header, $sequence);
    }

    $/ = $orig_record_sep; #reset the record separator to original setting.

    return ($seqobj); #returns null if not instantiated.
}


#### finish() closes the open filehandle to the query database.
sub finish {
    my $self = shift;
    my $filehandle = $self->{fileHandle};
    close $filehandle;
    $self->{fileHandle} = undef;
}

####
sub retrieve_all_seqs_hash {
	my $self = shift;
  my $cov= shift;
	my $acc_to_seq=();

	while (my $seq_obj = $self->next()) {
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();
    #next if($seq_obj->header->{KC})
		$acc_to_seq->{$acc} = $seq_obj;
	}

	return $acc_to_seq;
}

####
sub retrieve_all_seqs_array{
	my $self = shift;
	my $aseq=();
	my $index=0;
	while (my $seq_obj = $self->next()) {
		$seq_obj->{index}=$index;
		push(@{$aseq},$seq_obj);
		$index++;
	}
	return $aseq;
}

##############################################
package Sequence;
use strict;

sub new {
    my ($packagename, $header, $sequence) = @_;

    ## extract an accession from the header:
    my ($acc, $rest) = split (/\s+/, $header, 2);
    $sequence=uc($sequence);
    $sequence=~s/U/T/g;

    my $self = { accession => $acc,
		 header => $header,
		 sequence => $sequence,
     rsequence=> undef,
		 filename => undef,
         edges => undef, };
    #temporary to check the orientation in the bruijn grpah

    _bcalm_header($header,$self);
    _rev_comp($self,$sequence);
    bless ($self, $packagename);
    return ($self);
}

#aux funtion to parse the header of bcalm
sub _bcalm_header{
    my ($h,$hash)=@_;
    if($h=~/(\d+) LN:i:(\d+) KC:i:(\d+) km:f:(\d+)/){
        $hash->{LN}=$2;
        $hash->{KC}=$3;
        $hash->{KM}=$4;
        #$hash->{KCOV}=$3/$2;
        #$hash->{KMCOV}=$4/$2;
    }

    my @m = $h =~/L:([+,-]):(\d+):([+,-])/g;
    for(my $i=0;$i <= scalar(@m)-3; $i+=3 ){
            my $tmp=();
            $tmp->{o1}=$m[$i];
            $tmp->{uid}=$m[$i+1];
            $tmp->{o2}=$m[$i+2];
            #print join(" ",$tmp->{o1},$tmp->{uid},$tmp->{o2})."\n";
            push(@{$hash->{edges}},$tmp);
    }
}

sub _rev_comp{
  my ($h,$str)=@_;
  $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
   $h->{rsequence} = reverse $str;
}
####
sub get_accession {
    my $self = shift;
    return ($self->{accession});
}

###
sub get_length{
    my $self = shift;
    return length($self->{sequence});
}

####
sub get_header {
    my $self = shift;
    return ($self->{header});
}

sub set_header{
    my ($self,$header)=@_;
    $self->{header}=$header;

}
#return the reverse complement of current sequence
sub get_revcomp{
   my $self = shift;
   my $str=$self->{sequence};
   $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
   my $revseq = reverse $str;
   return $revseq;
}
####
sub get_sequence {
    my $self = shift;
    return ($self->{sequence});
}

####
sub get_FASTA_format {
    my $self = shift;
    my %settings = @_;

    my $fasta_line_len = $settings{fasta_line_len} || 60;

    my $header = join(" ",$self->get_accession(),$self->get_header());
    my $sequence = $self->get_sequence();
    if ($fasta_line_len > 0) {
        $sequence =~ s/(\S{$fasta_line_len})/$1\n/g;
        chomp $sequence;
    }
    my $fasta_entry = ">$header\n$sequence\n";
    return ($fasta_entry);
}


####
sub write_fasta_file {
    my $self = shift;
    my $filename = shift;

    my ($accession, $header, $sequence) = ($self->{accession}, $self->{header}, $self->{sequence});

	my $fasta_entry = $self->get_FASTA_format();

    my $tempfile;
    if ($filename) {
		$tempfile = $filename;
    } else {
		my $acc = $accession;
		$acc =~ s/\W/_/g;
		$tempfile = "$acc.fasta";
    }

    open (TMP, ">$tempfile") or die "ERROR! Couldn't write a temporary file in current directory.\n";
    print TMP $fasta_entry;
    close TMP;
    return ($tempfile);
}

####
sub get_core_read_name {
    my $self = shift;

    my $acc = $self->get_accession();
    $acc =~ s|/[12]$||;
    return($acc);
}



1; #EOM
