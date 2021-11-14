
###############################################################################
# Author: Carol Moraga Quinteros 
# Laboratory: ERABLE-UCBL
# Copyright (c)
# year: 2020
###############################################################################

use Data::Dumper;
use Getopt::Std;
use strict;

#create the unipath graph from the BCALM output
package UGraph;

use strict;
#fasta class that parse the BCALM output
use Fasta_reader;
#use Graph;
use Graph::Directed;
use Graph::Undirected;
use Data::Dumper;
#class that allow to compute some seq related stats
use SeqComplex;
sub new{
  my ($packagename, $bcalm_fasta, $cov) = @_;
  my $freader = new Fasta_reader($bcalm_fasta);
  my $unitigs=$freader->retrieve_all_seqs_hash();
  #print Dumper($unitigs);
  my $self = { unitigs => $unitigs,cov => $cov, freader=>$freader};

  bless ($self, $packagename);
  #method to rescue singletons miRNAs as potential candidates
  $self->_check_singletons();

  return ($self);
}



sub _check_singletons{
	my $self=shift;
	my $potential=();
	my $otherseqs=();
	my @methods = qw/gc gcs cpg cwf ce cm1 cm2 cm3 cm4 cm5 cm6 ct1 ct2 ct3 ct4 ct5 ct6 cl1 cl2 cl3 cl4 cl5 cl6/;
	foreach my $uni (sort { $self->{unitigs}->{$b}->{KC} <=> $self->{unitigs}->{$a}->{KC}} keys %{$self->{unitigs}}){
	    #we set the pointer to the Fasta_Reader object
   	 my $suni=$self->{unitigs}->{$uni};
	if( $suni->{edges} == undef and $suni->{KC} >=100 and $suni->{LN} <=23 and $suni->{LN} >=21 ){
	my $tmp=();
        $tmp->{id}=$suni->{accession};
        $tmp->{seq}=$suni->{sequence};
        $tmp->{LN}=$suni->{LN};
        $tmp->{KM}=$suni->{KM};
        $tmp->{KC}=$suni->{KC};
       	my %results = runAllMethods($tmp->{seq});
        #print Dumper(%results);
        #print O "$start_node";

        foreach my $m (@methods) {
      		my $val = shift @{ $results{$m} };
                $tmp->{comp2}->{$m}=$val;
      		print O "\t$val";
      	}
        #print Dumper($tmp);
        #we save the potential list of miRNAs to compare them
        push(@{$potential},$tmp);
	}
	#we save long singletons with high KM value as putative long non coding
	if( $suni->{edges} == undef and $suni->{KM} >=50 and $suni->{LN} >25 ){
        my $tmp=();
        $tmp->{id}=$suni->{accession};
        $tmp->{seq}=$suni->{sequence};
        $tmp->{LN}=$suni->{LN};
        $tmp->{KM}=$suni->{KM};
        $tmp->{KC}=$suni->{KC};
        my %results = runAllMethods($tmp->{seq});
        #print Dumper(%results);
        #print O "$start_node";

        foreach my $m (@methods) {
                my $val = shift @{ $results{$m} };
                $tmp->{comp2}->{$m}=$val;
                print O "\t$val";
        }
	#we save putative long non coding
	push(@{$otherseqs},$tmp);
	}
    }
    #print scalar(@{$potential})."\n";
    #print Dumper($potential);
	$self->{candidates_singletons}=$potential;
	$self->{candidates_otherseqs}=$otherseqs;

}


#we return a hash with the putative seed ids and the corresponding CC in the undirected graph
sub select_seeds{
  my $self = shift;
  #maximal seed length, should be set to the max miRNA length
  my $m_len= shift;
  $self->create_directed_graph();
  #min coverage for being selected as a potencial seed
  my $min_s_cov= shift;
  my $total_seeds=0;
  my $seeds=();
  foreach my $uni (sort { $self->{unitigs}->{$b}->{KC} <=> $self->{unitigs}->{$a}->{KC}} keys %{$self->{unitigs}}){
    #we set the pointer to the Fasta_Reader object
    my $suni=$self->{unitigs}->{$uni};
    #we skypt deleted unipaths
    next if(defined $self->{dunipaths}->{$uni});
    my $tmp=();
    #we remove seed that the length is larger than 50 and larger than 18
    if($suni->{LN} <= $m_len and $suni->{KM} >= $min_s_cov and $self->{sccc}->{$self->{scc}->{$uni}} < 10 and $self->{sccc}->{$self->{scc}->{$uni}} > 1){
    #print join(" ","SEED",$suni->get_accession(),$suni->{KC},$suni->{KM},$suni->{LN},$self->{scc}->{$uni},$self->{sccc}->{$self->{scc}->{$uni}})."\n";
    $total_seeds++;
    push(@{$seeds->{$self->{scc}->{$uni}}},$suni->get_accession());
    #print join(" ","SCCS",$suni->get_accession(),$self->{scc}->{$uni},$suni->{KC},$suni->{KM},$suni->{LN})."\n";
  }else{
    #we delete the unipath due that is not part of a SCC or live in a larger part of the graph
      $self->{dunipaths}->{$uni}++;
  }
  }
  #print Dumper($seeds);
  my $total_scc=scalar(keys %$seeds);
  print "Total potencial seeds $total_seeds  (Length <= $m_len and Coverage >= $min_s_cov) found in $total_scc Strongly Connected Components\n";
  #we save the seeds in the main graph structure
  $self->{seeds}=$seeds;
  $self->update_edges();
}

sub del_neigboords_by_coverage2{
  my $self=shift;
  my $max_cov_diff=shift;
  # from each unipath compares the coverage upto $depht_neighbourds ahead [1-3]
  my $depht_neighbourds=shift;

  my $del_edges=0;
  my $total_edges=0;
  my $del_unipath=();

  #we check the neigbourds of the same SCC
  foreach my $uni (sort { $self->{unitigs}->{$b}->{KC} <=> $self->{unitigs}->{$a}->{KC}} keys %{$self->{unitigs}}){
      #print join(" ",$s,$cc)."\n";
      my $suni=$self->{unitigs}->{$uni};
      #we add the edges
      #print join(" ",$cc,$s,$suni->{LN},$suni->{KM},$suni->{sequence})."\n";
      if(defined $suni->{edges}){
            my $tedges=();
            #level 1 of iterations
              foreach my $e(@{$suni->{edges}}){
                $total_edges++;
                my $nd1=$self->{unitigs}->{$e->{uid}};
                my $dc=0;
                my $break=0;
		next if (!defined $nd1);
                if($nd1->{KM} < $suni->{KM}){
                  $dc=$suni->{KM}/$nd1->{KM};
                  if($dc>= $max_cov_diff){
                  $del_unipath->{$nd1->get_accession()}++;
                   }
                }else{
                  $dc=$nd1->{KM}/$suni->{KM};
                  if($dc>= $max_cov_diff){
                   $del_unipath->{$suni->get_accession()}++;
                   }
                }

                if($dc>= $max_cov_diff){
                  $break=1;
                  $del_edges++;
                }else{
                  push(@{$tedges},$e);
                }
              }
              $suni->{edges}=undef;
              $suni->{edges}=$tedges;
       }

  }
  
#we check the deleted unipaths and as for number of remaining edges
  foreach my $uni(keys %{$del_unipath}){
    my $suni=$self->{unitigs}->{$uni};
    if(defined $suni->{edges}){
        #print join(" ","U",$suni->get_accession,$del_unipath->{$uni},scalar(@{$suni->{edges}}))."\n";
        #the unipath as conecctions to others unipahts
        my $deletedu=0;
        my $tedges=0;
        foreach my $e(@{$suni->{edges}}){
          if(defined $del_unipath->{$e->{uid}}){
            $deletedu++;
          }
          $tedges++;
        }
        if($deletedu == $tedges){
            $self->{dunipaths}->{$uni}++;
        }

    }else{
      $self->{dunipaths}->{$uni}++;
    }
  }
  print "A total of $del_edges edges ($total_edges) were removed due to huge coverage differences [>=$max_cov_diff] $del_edges/$total_edges\n";
  $self->update_edges();
}

#we remove some topologies related to sequencing errors
sub delete_by_topology{
      my $self= shift;
      my $prefix= shift;
      open(DELETE,">".$prefix.".del_by_topology.txt") or die "cannot create topology file\n";
      #we create the directed graph and compute SCC and SCC node count
      $self->create_directed_graph();
      my $cases=0;
      my $cases3=0;
      my $cases4=0;
      my $cases5=0;
      my $putatives=0;
      foreach my $scc (keys %{$self->{sccc}}){
        #putatives miRNAs
        if($self->{sccc}->{$scc} <= 10 and $self->{sccc}->{$scc} >= 2){
          $putatives++;
        }
        #we work with SCC with total number of unipaths of 3,4,5 and we attempt to remove them
        next if($self->{sccc}->{$scc} <= 2 or $self->{sccc}->{$scc} > 5);
        #print join(" ",$scc,$self->{sccc}->{$scc})."\n";
        #print Dumper($self->{scc2nodes}->{$scc});
        #We iterate the nodes of the SCC
        my $degree=();
        foreach my $n (@{$self->{scc2nodes}->{$scc}}){
                my $uni=$self->{unitigs}->{$n};
                #print Dumper($uni);
                #print join(" ",$uni->get_accession,$uni->{KC},$uni->{LN},$uni->{KM})."\n";
                my $inout=();
                foreach my $e(@{$uni->{edges}}){
                         #print join(" ",$uni->get_accession,$e->{o1},$e->{uid},$e->{o2})."\n";
                         #print join(" ",$uni->get_accession,$e->{o1},$scc)."\n";
                         $degree->{$uni->get_accession}->{$e->{o1}}++;
                }
                if(!defined $degree->{$uni->get_accession}->{"+"}){
                  $degree->{$uni->get_accession}->{"+"}=0;
                }
                if(!defined $degree->{$uni->get_accession}->{"-"}){
                  $degree->{$uni->get_accession}->{"-"}=0;
                }
        }
        #print "\n";
        my $deleteG=0;
        #case 3
        if(scalar(keys %{$degree}) == 3){
          foreach my $u (keys %{$degree}){
              #print Dumper($degree->{$u});
                 foreach my $ori (keys %{$degree->{$u}}){
                    if($degree->{$u}->{$ori} == 2){
                      $deleteG=1;
                      #$degree->{delete3}=1;
                      $cases3++;
                    }
                 }
           }
        }
        #case 4
        if(scalar(keys %{$degree}) == 4){
          my $counter4=0;
          my $counter4_1=0;
          foreach my $u (keys %{$degree}){
                 foreach my $ori (keys %{$degree->{$u}}){
                    if($degree->{$u}->{$ori} == 2){
                      $counter4++;
                    }
                    if($degree->{$u}->{$ori} == 3){
                      $counter4_1=1;
                    }
                }
           }
           #we delete if at least two unipaths as the same in or out degree
           if($counter4 >= 2 or $counter4_1 == 1){
             $deleteG=1;
             $cases4++;
           }
        }
        #case 5
        if(scalar(keys %{$degree}) == 5){
          my $counter5=0;
          my $counter5_3=0;
          my $counter5_4=0;
          my $counter5_2=0;
          my $counter5_5=0;
          foreach my $u (keys %{$degree}){

                 foreach my $ori (keys %{$degree->{$u}}){
                    if($degree->{$u}->{$ori} == 2){
                      $counter5++;
                    }


                    if($degree->{$u}->{$ori} == 3){
                      $counter5_5++;
                    }
                }

                if(($degree->{$u}->{"+"}==2 and $degree->{$u}->{"-"} ==1) or ($degree->{$u}->{"+"}==1 and $degree->{$u}->{"-"} ==2)){
                    $counter5_3=1;
                }
                #we ask for 2 and 2
                if($degree->{$u}->{"+"}==2 and $degree->{$u}->{"-"} ==2){
                  $counter5_4=1;
                }
                if($degree->{$u}->{"+"}==1 and $degree->{$u}->{"-"} ==1){
                  $counter5_2=1;
                }
           }
           #we delete if at least two unipaths as the same in or out degree
           #if($counter5 == 5 or ( $counter5=1 and $counter5_3 == 1) or ($counter5_4 == 1)){
           if($counter5_4 == 1 or $counter5_3 == 1 or ($counter5 ==2 and $counter5_2 == 1) or ($counter5==1 and $counter5_5==1) or $counter5_5 > 1){
             $deleteG=1;
             $cases5++;
           }
        }

        #delete the unipaths from the graph
        if($deleteG){
            foreach my $u (keys %{$degree}){
              $self->{dunipaths}->{$u}=1;
	      print DELETE $u."\n";
              #we delete the reference of the node to the SCC
              delete $self->{scc}->{$u};
            }
            #we delete the scc2nodes
            delete $self->{scc2nodes}->{$scc};
            #we delete the sccc counter (number of unipaths in SCC)
            delete $self->{sccc}->{$scc};
            #we delete the
        }

        $cases++;
      }
      my $total_deleted_topology=$cases3+$cases4+$cases5;
      print "Potential $putatives delete by topology3(Y): $cases3 (".($cases3/$putatives * 100)."%) \n";
      print "Potential $putatives delete by topology4(Y): $cases4 (".($cases4/$putatives * 100)."%) \n";
      print "Potential $putatives delete by topology5(Y): $cases5 (".($cases5/$putatives * 100)."%) \n";
      print "Potential $putatives delete by topology: $total_deleted_topology (".($total_deleted_topology/$putatives * 100)."%) \n";
      $self->update_edges();
      close(DELETE);
      #print Dumper($self->{scc2nodes});
}

#we assemble each SCC into its unipaths and classify each SCC into protential miRNA or other kind of sequences
sub assemble_unipaths{
    #UGRAPH object
    my $self=shift;
    #kmer_size
    my $kmer_size=shift;
    #the lenght of edges
    $kmer_size--;
    my $prefix=shift;
    #seq related methods
    my @methods = qw/gc gcs cpg cwf ce cm1 cm2 cm3 cm4 cm5 cm6 ct1 ct2 ct3 ct4 ct5 ct6 cl1 cl2 cl3 cl4 cl5 cl6/;
    open O, ">$prefix.seq_composition_and_complexity.txt" or die "cannot open $prefix.seq_composition_and_complexity.txt\n";
    open LARGESEQ, ">$prefix.other_sequences.txt" or die "cannot open $prefix.other_sequences.txt\n";
    open ASM, ">$prefix.assembled_sequences.txt" or die "cannot open $prefix.assembled_sequences.txt\n";
    open CC, ">$prefix.discarded_ccs.txt" or die "cannot open $prefix.discarded_ccs.txt\n";
    open BRA, ">$prefix.discarded_branching.txt" or die "cannot open $prefix.discarded_branching.txt\n";

    # Print header
    print O "seq";
    foreach my $m (@methods) { print O "\t$m"; }
    print O "\n";
    my $potencial=();
    my $otherseqs=();
    my $putatives=0;
    foreach my $scc (keys %{$self->{sccc}}){
      #putatives miRNAs
      if($self->{sccc}->{$scc} <= 10 and $self->{sccc}->{$scc} >= 2){
        $putatives++;
      }else{
        print CC join(" ",$scc, join(",",@{$self->{scc2nodes}->{$scc}}))."\n";
      }
      my $degree=();
      #we iterate the nodes of the SCC and we compute the indegree/outdegree to determine branching nodes and head/tail nodes
      my $head=();#the node
      my $branching=();
      my $internal=0;
      my $head_tail=0;
      my $branching_count=0;
      my $scc_km=0;
      my $scc_kc=0;
	my $uni=();
	my $cont=0;
      foreach my $n (@{$self->{scc2nodes}->{$scc}}){
               $uni=$self->{unitigs}->{$n};
              foreach my $e(@{$uni->{edges}}){
                       $degree->{$uni->get_accession}->{$e->{o1}}++;
              }
              if(!defined $degree->{$uni->get_accession}->{"+"}){
                $degree->{$uni->get_accession}->{"+"}=0;
              }
              if(!defined $degree->{$uni->get_accession}->{"-"}){
                $degree->{$uni->get_accession}->{"-"}=0;
              }
              #we determine the number of internal nodes, which are the ones that shold be reduced
              if($degree->{$uni->get_accession()}->{"+"} == 1 and $degree->{$uni->get_accession}->{"-"} == 1){
                $internal++;
              #head/tail node
              }elsif(($degree->{$uni->get_accession()}->{"+"} == 0 and $degree->{$uni->get_accession}->{"-"} == 1) or ($degree->{$uni->get_accession()}->{"+"} == 1 and $degree->{$uni->get_accession}->{"-"} == 0) ){
                $head_tail++;
                $head->{$uni->get_accession()}=1;
              }else{
                   #is a branching node, one that stop the unipath extention
                   $branching_count++;
                   $branching->{$uni->get_accession}=1;
              }
              $scc_km+=$uni->{KM};
              $scc_kc+=$uni->{KC};
       }
       #we save singleton nodes with km value more than 20 ---- AQUI
       if(scalar(keys %{$degree})==1){
	#my $uni=$self->{unitigs};

        if( $uni->{KM} >=50 and $uni->{LN} <=25  and $uni->{LN} >=18){
        	my $tmp=();
        	$tmp->{id}=$uni->{accession};
        	$tmp->{seq}=$uni->{sequence};
        	$tmp->{LN}=$uni->{LN};
        	$tmp->{KM}=$uni->{KM};
        	$tmp->{KC}=$uni->{KC};
        	my %results = runAllMethods($tmp->{seq});
		$cont++;

        foreach my $m (@methods) {
                my $val = shift @{ $results{$m} };
                $tmp->{comp2}->{$m}=$val;
                print O "\t$val";
        }
	print O "\n";
	#print Dumper($uni);
        #print Dumper($tmp);
    	#print Dumper($potential)
        #we save the potential list of miRNAs to compare them
        push(@{$self->{candidates_singletons}},$tmp);
        }
     }
       if($branching_count > 0 or $head_tail==0){

         foreach my $n (@{$self->{scc2nodes}->{$scc}}){
                  $uni=$self->{unitigs}->{$n};
                  print BRA join(" ",$scc,$uni->get_accession,$uni->{KM},$uni->{KC},$uni->{sequence},join(",",@{$uni->{edges}}))."\n";

         }
         next;
       }



       # we compact the unipath
       my ($start_node,$end_node)=sort keys %{$head};

       #print join(" ",  $scc,$start_node,$end_node)."\n";
       my $unodes=();
       my $uexp=();
       my $e=$self->{unitigs}->{$start_node}->{edges}[0];
       #print $self->{unitigs}->{$start_node}->{sequence}."\n";
       my $sequni="";
       push(@{$unodes},$start_node);
       push(@{$uexp},$self->{unitigs}->{$start_node}->{KM});

       while($start_node != $end_node){
         my $prev=$self->{unitigs}->{$start_node};
          my $current=$self->{unitigs}->{$e->{uid}};
          push(@{$unodes},$current->get_accession);
          push(@{$uexp},$current->{KM});
          #unipath is empty
          if(length($sequni) eq 0){
            $sequni = ($e->{o1} eq "+") ? $prev->{sequence}:$prev->{rsequence};
          }
          $sequni .= ($e->{o2} eq "+") ? substr($current->{sequence},$kmer_size,$current->{LN}-$kmer_size):substr($current->{rsequence},$kmer_size,$current->{LN}-$kmer_size);
          #we current update the edge
          $e=$current->{edges}[0];
          if($e->{uid} == $start_node){
              $e=$current->{edges}[1];
          }
          $start_node=$current->get_accession;

       }
       #print $sequni."\n";
       print ASM "Assembled_unipaths_nodes:".join(",",@{$unodes})." KM:".join(",",@{$uexp})."SEQ:$sequni L:".length($sequni)."\n";
       #if(length($sequni) > 25){
       #  check_cov_asm($sequni,$unodes,$uexp);
       #}
       #check if we have to trim
       if(length($sequni) >=18 and length($sequni)<=25){
        #print "Potencial miRNA: $start_node SEQ=$sequni L=".length($sequni)." KM=".$scc_kc/(length($sequni)-$kmer_size)."\n";
        my $tmp=();
        $tmp->{id}=$start_node;
        $tmp->{seq}=$sequni;
        $tmp->{LN}=length($sequni);
        $tmp->{KM}=$scc_kc/(length($sequni)-$kmer_size);
        $tmp->{KC}=$scc_kc;
       	my %results = runAllMethods($sequni);
        #print Dumper(%results);
        print O "$start_node";

        foreach my $m (@methods) {
      		my $val = shift @{ $results{$m} };
          $tmp->{comp2}->{$m}=$val;
      		print O "\t$val";
      	}
        #print Dumper($tmp);
        #we save the potencial list of miRNAs to compare them
        push(@{$potencial},$tmp);
        print O "\n";
        #unipaths larger than 50 are saved in the LARGESEQ file
      }elsif(length($sequni)>=26){
        print LARGESEQ ">$start_node L=".length($sequni)." COVERAGE=".$scc_kc/(length($sequni)-$kmer_size)."\n";
        print LARGESEQ $sequni."\n";
        my $tmp=();
        $tmp->{id}=$start_node;
        $tmp->{seq}=$sequni;
        $tmp->{LN}=length($sequni);
        $tmp->{KM}=$scc_kc/(length($sequni)-$kmer_size);
        $tmp->{KC}=$scc_kc;
        push(@{$otherseqs},$tmp);
      }

    }
	#add singletons as potential candidates
      foreach my $us(@{$self->{candidates_singletons}}){
	 push(@{$potencial},$us);
      }
      #save the potencial mirRNAS
      $self->{potential_miRNA}=$potencial;
       #we print otherseqs coming from  unipaths
      foreach my $ouni (@{$self->{candidates_otherseqs}}){
		print LARGESEQ ">$ouni->{id} L=$ouni->{LN} COVERAGE=$ouni->{KM}\n";
		print LARGESEQ $ouni->{seq}."\n";
		push(@{$otherseqs},$ouni);
      }
      #save the assembled and longer sequences
      $self->{otherseqs}=$otherseqs;
      #print Dumper($self->{potential_miRNA});
      $self->update_edges();
      close(O);
      close(LARGESEQ);
      close (ASM);
      close(BRA);
      close(CC);
}


sub check_cov_asm{
    my ($seq,$unis,$cov)=@_;
    #check the max KM
    my $max=0;
    my $max_pos=0;

    for(my $i=0; $i<@{$cov}; $i++){
          if(@$cov[$i] > $max){
            $max=@$cov[$i];
            $max_pos=$i;
          }
    }

    my $split=0;
    for (my $i=0;$i<@{$cov}-1;$i++){
         if(@$cov[$i] < $max/2){
           $split++;
           print join(" ","D",@$cov[$i],$max/2,$i,@$unis[$i])."\n";
         }
    }
    if($split){
      print "\n\n";
      #we iterate the sequence from start to end
      print join(" ",length($seq),$seq,$max)."\n";
      print join(" ",@{$cov})."\n";
      print join(" ",@{$unis})."\n";
      print "\n\n";
    }


}


#it the edit distance to compute overlap among the long reads
sub compute_overlap_candidates_edit{
   my $self= shift;
   #length of the overlap size
   my $o_size = shift;
   #file to save the cluster in fasta format
   my $prefix = shift;
   #path to ed_aligner
   my $aln = shift;

   my $file=">".$prefix.".edit_clusters.ED".$o_size.".txt";
   #remove the overlap
    open(CLUSTER,$file) or die "cannot open file $file\n";
   #we create an undirected graph
   my $g2 = Graph::Undirected->new; #we create an undirected graph
   my $kmer_overlaps=();
   print join(" ","Potencial miRNA before edit distance  =" , scalar(@{$self->{potential_miRNA}}))."\n";
   my $hash2pos=();
   my $i=0;
   my $fatmp=">".$prefix.".ed.fa";

    open(FATMP,$fatmp) or die "cannot open file $fatmp\n";
  foreach my $c(@{$self->{potential_miRNA}}){
        #print join(" ",$c->{id},$c->{seq})."\n";
        print FATMP join("\n",">".$c->{id},$c->{seq})."\n";
        $g2->add_vertex($c->{id});
        $hash2pos->{$c->{id}}=$i;
        $i++;

     }
  close(FATMP);
  # we call ed_aliner to compute the overlaps among the candidates
   my $cmd="$aln -m SHW -k $o_size -p -f NICE $prefix.ed.fa $prefix.ed.fa > $prefix.ed.log";
   #my $cmd="$aln -m NW -k $o_size -p -f NICE $prefix.ed.fa $prefix.ed.fa > $prefix.ed.log";
   print $cmd."\n";
  system($cmd) == 0 or die "cannot run ed_aligner";
  unlink $prefix.".ed.fa";
  open(EDR,$prefix.".ed.log") or die "cannot open file $prefix.ed.log\n";

  my $miRNAS=();
  while(my $line=<EDR>){
        chomp $line;
	next if($line !~/^#/);
	$line=~s/#//;
	$line=~s/:/ /;
        my @data=split / /,$line;
	#we skypt equal candidates
	next if($data[0] == $data[1]);
        #print Dumper(@data);
	#print join(" ",$data[0],$data[1])."\n";
	#printf("#%s %s: %d %d %d %d %d\n", ntarget, queries[i].name.c_str(), scores[i], numLocations[i],targetLength,queryLength,alignmentLength);
	#we add the edges to the graph
        $g2->add_edge($data[0],$data[1]);
	$miRNAS->{$data[0]}++;
	$miRNAS->{$data[1]}++;
  }
  close(EDR);

     my @cc=$g2->connected_components();
     print join(" ","Total miRNAs whit overlaps  at edit distance (ED=$o_size)",scalar(keys %{$miRNAS}))."\n";
     print "Total CC=".scalar(@cc)."\n";
     my $cc_index=0;
     my $pass_filter=();
     foreach my $c(@cc){
       my $selected=0;
       my $max=0;
       foreach my $n(@{$c}){
         my $node=$self->{potential_miRNA}[$hash2pos->{$n}];
         if($node->{KM} > $max){
           $selected=$node->{id};
           $max=$node->{KM};
         }
         #print Dumper($node);
         #'KM' => '304',
        #  'id' => 70984,
        #'LN' => 20,
        #'seq' => 'ACTTTTTCTCTGACCAACCA',
         print CLUSTER join(" ",$cc_index,$n,$node->{id},$node->{KM},$node->{LN},$node->{seq})."\n";
       }
       #print join(" ","***",$selected,$max)."\n";
       #we delete the non selected nodes from the array
       my $snode=$self->{potential_miRNA}[$hash2pos->{$selected}];
       print CLUSTER join(" ","***",$cc_index,$selected,$snode->{id},$snode->{KM},$snode->{LN},$snode->{seq})."\n\n";
       push(@{$pass_filter},$snode);
       $cc_index++;
       #print "\n";
     }
     #we sort the array by pass  km
     $self->{potential_miRNA}=$pass_filter;
     print join(" ","Selected nodes",scalar(@{$self->{potential_miRNA}}))."\n";
     close(CLUSTER);
}




#function that compute overlap ammong the selected miRNA candidates
#it use a hash to compute the overlap among them
sub compute_overlap_candidates_exact{
   my $self= shift;
   #length of the overlap size
   my $o_size = shift;
   #file to save the cluster in fasta format
   my $prefix = shift;
   my $file=">".$prefix.".exact_clusters.k".$o_size.".txt";
   #remove the overlap
    open(CLUSTER,$file) or die "cannot open file $file\n";
   #we create an undirected graph
   my $g2 = Graph::Undirected->new; #we create an undirected graph
   my $kmer_overlaps=();
   if(!defined $self->{potential_miRNA}){
		print "No candidate founds\n";
		exit(0);
   }
   print join(" ","Potencial miRNA before overlaps  =" , scalar(@{$self->{potential_miRNA}}))."\n";
   my $hash2pos=();
   my $i=0;
  foreach my $c(@{$self->{potential_miRNA}}){
        #print join(" ",$c->{id},$c->{seq})."\n";
        $g2->add_vertex($c->{id});
        $hash2pos->{$c->{id}}=$i;
        $i++;
        #hashing in foward strand
        my $numKmersInRead=length($c->{seq})-$o_size;
        #hashing in reverse strand
        my $revseq=$c->{seq};
        $revseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
        for(my $j=0;$j<$numKmersInRead;$j++){
           $kmer_overlaps->{substr($c->{seq},$j,$o_size)}->{$c->{id}}++;
           $kmer_overlaps->{substr($revseq,$j,$o_size)}->{$c->{id}}++;
        }
     }

     my $miRNAS=();
     foreach my $k (keys %{$kmer_overlaps}){
       next if(scalar(keys %{$kmer_overlaps->{$k}}) == 1);
       #print join(" ","OVERLAPS",$k, scalar(keys %{$kmer_overlaps->{$k}}))."\n";
       my @nodes=keys %{$kmer_overlaps->{$k}};
       $miRNAS->{$_}++ foreach(@nodes);
       for(my $i=0;$i <= $#nodes; $i++){
         for(my $j=$i+1; $j <=$#nodes; $j++){
           next if($nodes[$i] == $nodes[$j]);
            $g2->add_edge($nodes[$i],$nodes[$j]);
           #print join(" ","Addging edge between",$nodes[$i],$nodes[$j])."\n";
         }
       }
     }
     my @cc=$g2->connected_components();
     print join(" ","Total miRNAs whit overlaps (L=$o_size)",scalar(keys %{$miRNAS}))."\n";
     print "Total CC=".scalar(@cc)."\n";
     my $cc_index=0;
     my $pass_filter=();
     foreach my $c(@cc){
       my $selected=0;
       my $max=0;
       foreach my $n(@{$c}){
         my $node=$self->{potential_miRNA}[$hash2pos->{$n}];
         if($node->{KM} > $max){
           $selected=$node->{id};
           $max=$node->{KM};
         }
         print CLUSTER join(" ",$cc_index,$n,$node->{id},$node->{KM},$node->{LN},$node->{seq})."\n";
       }
       #we delete the non selected nodes from the array
       my $snode=$self->{potential_miRNA}[$hash2pos->{$selected}];
       print CLUSTER join(" ","***",$cc_index,$selected,$snode->{id},$snode->{KM},$snode->{LN},$snode->{seq})."\n\n";
       push(@{$pass_filter},$snode);
       $cc_index++;
       #print "\n";
     }
     #we sort the array by pass  km
     $self->{potential_miRNA}=$pass_filter;
     print join(" ","Selected nodes",scalar(@{$self->{potential_miRNA}}))."\n";
     close(CLUSTER);
}

#funtion that overlap candidate miRNA to the longer detected sequences
sub overlap_to_othersequences{
  my $self=shift;

  my $o_size=shift;
  my $prefix=shift;

  open (OVL,">".$prefix.".ovl_to_long_seqs.txt") or die "cannot create file\n";
  #print Dumper($self->{otherseqs});
  my $kmer_overlaps=();
  print join(" ","Potencial miRNA before overlaps with long assembled sequences  =" , scalar(@{$self->{potential_miRNA}}))."\n";
  my $hash2pos=();
  my $i=0;
 foreach my $c(@{$self->{otherseqs}}){
       #print join(" ",$c->{id},$c->{seq})."\n";
       #$g2->add_vertex($c->{id});
       $hash2pos->{$c->{id}}=$i;
       $i++;
       #hashing in foward strand
       my $numKmersInRead=length($c->{seq})-$o_size;
       #hashing in reverse strand
       my $revseq=$c->{seq};
       $revseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
       for(my $j=0;$j<$numKmersInRead;$j++){
          $kmer_overlaps->{substr($c->{seq},$j,$o_size)}->{$c->{id}}++;
          $kmer_overlaps->{substr($revseq,$j,$o_size)}->{$c->{id}}++;
       }
    }
    #print Dumper($kmer_overlaps);
    my @scandidates = sort {$b->{KM} <=> $a->{KM}} @{$self->{potential_miRNA}};
    my $ii=0;
    my $pass_filter=();
    foreach my $mi(@scandidates){
      #print join(" ","MIP",$mi->{id},int($mi->{KM}),$mi->{LN},$mi->{seq},$ii)."\n";
      #hashing in foward strand
      my $numKmersInRead=length($mi->{seq})-$o_size;
      #hashing in reverse strand
      my $revseq=$mi->{seq};
      $revseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
      my $overlap=0;
      for(my $j=0;$j<$numKmersInRead;$j++){
         if(defined $kmer_overlaps->{substr($mi->{seq},$j,$o_size)}){
            #print join(" ","OVERLAP",$mi->{id},substr($mi->{seq},$j,$o_size), join(",",keys %{$kmer_overlaps->{substr($mi->{seq},$j,$o_size)}}))."\n";
            $overlap=1;
         }
      }
      $ii++;
      #do not overlap long assembled sequences
      if($overlap == 0){
        push(@{$pass_filter},$mi);
      }else{
          print OVL join(" ","OVL_TO_LONG_SEQ($o_size):",$mi->{id},$mi->{seq},$mi->{KM},$mi->{LEN})."\n";
          push(@{$pass_filter},$mi);
      }
    }

    #we sort the array by pass  km
    $self->{potential_miRNA}=$pass_filter;
    print join(" ","Selected nodes",scalar(@{$self->{potential_miRNA}}))."\n";
    close (OVL);
}


#function that match kmers found in the RFAM database to the candidates miRNA
sub overlap_candidates_rfam{
  my $self=shift;
  my $o_size=shift;
  my $rfamfile=shift;#kmer database of rfam shold be the same size of k
  my $prefix=shift;
  my $kmer_overlaps=();
  print join(" ","Indexing potencial miRNAs for lookup in RFAM  =" , scalar(@{$self->{potential_miRNA}}))."\n";
  my $i=0;
 foreach my $c(@{$self->{potential_miRNA}}){
       $i++;
       #hashing in foward strand
       my $numKmersInRead=length($c->{seq})-$o_size;
       #hashing in reverse strand
       my $revseq=$c->{seq};
       $revseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
       for(my $j=0;$j<$numKmersInRead;$j++){
          $kmer_overlaps->{substr($c->{seq},$j,$o_size)}->{miRNA}->{$c->{id}}++;
          $kmer_overlaps->{substr($c->{seq},$j,$o_size)}->{count}=0;
          $kmer_overlaps->{substr($revseq,$j,$o_size)}->{miRNA}->{$c->{id}}++;
          $kmer_overlaps->{substr($revseq,$j,$o_size)}->{count}=0;
       }
  }
  #print Dumper($kmer_overlaps);
  open(RFAM,"gzip -dc $rfamfile |") or die "cannot open $rfamfile\n";
  my $match=0;
  while (my $line=<RFAM>) {
    # body...
    chomp $line;
    my ($kmer,$count)=split(" ",$line);
    #print join(" ",$kmer,$count)."\n";
    if(defined $kmer_overlaps->{$kmer}){
      $match++;
      $kmer_overlaps->{$kmer}->{count}+=$count;
    }

  }

  print "A total of $match kmers from RFAM as match with the candidate miRNAs (L=$o_size)\n";
  my $hits=();
  foreach my $k(keys %{$kmer_overlaps}){
    if($kmer_overlaps->{$k}->{count} > 0){
      #print join(" ",$k,$kmer_overlaps->{$k}->{count},join(",",keys %{$kmer_overlaps->{$k}->{miRNA}}))."\n";
      $hits->{$_}+=$kmer_overlaps->{$k}->{count} foreach(keys %{$kmer_overlaps->{$k}->{miRNA}});
    }
  }

  print join(" ","A total of ",scalar(keys %{$hits}), "miRNA candidates as hits with RFAM database")."\n";
  my $results=();
  open(RFAMH,">".$prefix.".RFAM_HITS.txt") or die "cannot create RFAM hits file\n";
  print RFAMH join(" ","#UID","KM","LN","SEQ","KHITS_RFAM")."\n";
  foreach my $c(@{$self->{potential_miRNA}}){
    if(defined $hits->{$c->{id}}){
      print RFAMH join(" ",$c->{id},$c->{KM},$c->{LN},$c->{seq},$hits->{$c->{id}})."\n";
    }else{
      push(@{$results},$c);
    }
  }

  close(RFAMH);
  $self->{potential_miRNA}=$results;
}

#output a fasta file and a table file with other features
sub write_potencial_miRNA{
  my $self=shift;
  my $prefix=shift;
  my $minkm=shift;

  open(FASTA,">".$prefix.".candidate_miRNA.fasta") or die "cannot create the candidate miRNA fasta file\n";
  open(TAB,">".$prefix.".candidate_miRNA.tsv") or die "cannot create the candidate miRNA tsv file\n";

  my @scandidates = sort {$b->{KM} <=> $a->{KM}} @{$self->{potential_miRNA}};
  my $rank=1;
  my @methods = qw/gc gcs cpg cwf ce/;
  print TAB join(" ","#UID","RANK","KM","KC","LN","GC","GCS","CPG","CWF","CE","SEQ")."\n";
  foreach my $c(@scandidates){
	next if($c->{KM} < $minkm);
      #we print a fasta file like
      print FASTA join(" ",">".$c->{id},"RANK=".$rank,"KM=".$c->{KM},"KC=".$c->{KC},"LN=".$c->{LN})."\n";
      print FASTA $c->{seq}."\n";
      #we print the table with the seq_composition_and_complexity values
      my @comp2=();
      push(@comp2,$c->{comp2}->{$_}) foreach(@methods);
      print TAB join(" ",$c->{id},$rank,$c->{KM},$c->{KC},$c->{LN},@comp2,$c->{seq})."\n";
      $rank++;
  }
 print "A total of $rank candidates were found\n";
}


sub create_directed_graph{
    my $self = shift;
    #we create a bidirected directed graph
    my $g1 = Graph::Directed->new;
    my $nodes=0;
    my $edges=0;
    #First we add al the vertices of the Graph
    foreach my $uni (keys %{$self->{unitigs}}){
       my $suni=$self->{unitigs}->{$uni};
       #we skypt deleted unipaths
       next if(defined $self->{dunipaths}->{$uni});
       #we add the node to the directed graph
       $g1->add_vertex($uni);
       $nodes++;
     }
     #we add the edges from the computation
     foreach my $uni (keys %{$self->{unitigs}}){
        my $suni=$self->{unitigs}->{$uni};
        #we skypt deleted unipaths
        next if(defined $self->{dunipaths}->{$uni});
        #the unipath dont have edges
        next if(!defined $suni->{edges});
        foreach my $e(@{$suni->{edges}}){
            #print join(" ","E",$e->{o1},$e->{uid},$e->{o2})."\n";
            #next if($e->{uid})
            #the node was deleted
            next if(defined $self->{dunipaths}->{$e->{uid}});
            $g1->add_edge($uni,$e->{uid});
            $edges++;
        }
     }
     #we save the graph
     $self->{bidirected}=$g1;
     #we compute the strongly connected components, where from each node is possible to reach any other
      my @scc = $g1->strongly_connected_components;
      my $index=0;
      my $hscc=();
      #the component length
      my $hsccc=();
      #the components to nodes
      my $scc2nodes=();
      foreach my $sc(@scc){
          foreach my $n (@{$sc}){
            #print join(" ",$index,$n)."\n";
            $hscc->{$n}=$index;
            $hsccc->{$index}++;
          }
          $scc2nodes->{$index}=$sc;
          $index++;
      }
      print "Total number of SCC $index\n";
     $self->{scc}=$hscc;
     $self->{sccc}=$hsccc;
     $self->{scc2nodes}=$scc2nodes;
     print "A total of $nodes nodes and A total of $edges edges were addeed to the bi-directed graph\n";
}


#remove  erroneus unipaths from  graph
#We search the nodes that do not have an incoming or outgoing edge
sub remove_tips{
    my $self= shift;
    #maximal seed length, should be set to the max miRNA length
    my $m_len= shift;
    #min coverage for being selected as a potencial seed
    my $min_s_cov= shift;
    my $removed_isolated_unipaths=0;
    my $tip_d1=0;
    my $total_unipaths=0;
    my $runi=();
    # we have to compute the incomming and outcomming degree
    foreach my $uni (sort { $self->{unitigs}->{$b}->{KC} <=> $self->{unitigs}->{$a}->{KC}} keys %{$self->{unitigs}}){
                my $suni=$self->{unitigs}->{$uni};
                $total_unipaths++;
                if(defined $suni->{edges}){
                  #is a tip, I have to compute the indegree and outdegree of each tip
                  if(scalar(@{$suni->{edges}}) == 1 and ($suni->{LN} > $m_len or $suni->{KM} < $min_s_cov)){
                  #print join(" ",$suni->get_accession(),$suni->{KC},$suni->{KM},$suni->{LN},scalar(@{$suni->{edges}}))."\n";
                   $tip_d1++;
                   #we save the unipath id in the list of removed unipahts
                   $runi->{$suni->get_accession()}++;
                 }else{
                    #we compute the indegree and out degree of the kmer
                    #means that the k-mer as conexions in both sides
                    next if(scalar(@{$suni->{edges}}) > 4);
                    #print join(" ",$suni->get_accession(),$suni->{KC},$suni->{KM},$suni->{LN},scalar(@{$suni->{edges}}))."\n";
                    my $inout=();
                    foreach my $e(@{$suni->{edges}}){
                        #print join(" ","E",$e->{o1},$e->{uid},$e->{o2})."\n";
                        $inout->{$e->{o1}}++;
                    }
                    #means that have unipaths coming in/out
                    next if(scalar(keys %{$inout}) > 1);
                    if($suni->{LN} > $m_len or $suni->{KM} < $min_s_cov){
                      #print join(" ",$suni->get_accession(),$suni->{KC},$suni->{KM},$suni->{LN},scalar(keys %{$inout}))."\n";
                      $runi->{$suni->get_accession()}++;
                      $tip_d1++;
                   }

                  }

                }else{
                  #isolated unipath, we check if the length and the coverage is correct
                  if($suni->{LN} > $m_len or $suni->{KM} < $min_s_cov){
                    #print join(" ","Removed isolated node",$suni->get_accession(),$suni->{KC},$suni->{KM},$suni->{LN},0)."\n";
                    $removed_isolated_unipaths++;
                    $runi->{$suni->get_accession()}++;
                  }
                }

      }
      #we have to remove the deleted nodes from as well as from the edges
      $self->{dunipaths}=$runi;
      #print join(" ")
      print "Total unipaths $total_unipaths\n";
      print "A total of $removed_isolated_unipaths (".$removed_isolated_unipaths/$total_unipaths * 100 .") isolated unipaths were removed from the graph (Length>$m_len or Coverage<$min_s_cov)\n";
      print "A total of $tip_d1 (".$tip_d1/$total_unipaths * 100 .") tips with degree 1 were removed from the graph (Length > $m_len or Coverage<$min_s_cov)\n";
      my $unipaths_alive= $total_unipaths - ($removed_isolated_unipaths+$tip_d1);
      print "A total of $unipaths_alive (". ($unipaths_alive/$total_unipaths) *   100 ." ) unipahts are alive\n";
      #print Dumper(scalar(keys %{$self->{dunipaths}}));
      $self->update_edges();
}





#methods that remove the removed edges from the graph
#useful for update the grap structure
sub update_edges{
    my $self = shift;
    # my @deleted=(46153,46204,109604,139159,46155,194324,252799,84563,119670,242338,281226,2509,172683,268972,272449);
    #    foreach my $u(@deleted){
		# my $suni=$self->{unitigs}->{$u};
		#if(defined $self->{dunipaths}->{$u}){
		#	print $u." was deleted\n";
		#}else{
		#	if(defined $suni){
		#	my @ee=();
		#	push(@ee,$_->{uid}) foreach(@{$suni->{edges}});
		#	 print join(" ",$u,$suni->{LN},$suni->{KM},$suni->{KC},"N_Edges=".scalar(@{$suni->{edges}}), join(",",@ee))."\n";
		#	}else{
		#	  print $u." is deleted\n";
		#	}
  #}
	#}
    foreach my $uni (sort { $self->{unitigs}->{$b}->{KC} <=> $self->{unitigs}->{$a}->{KC}} keys %{$self->{unitigs}}){
              my $suni=$self->{unitigs}->{$uni};
              if(defined $self->{dunipaths}->{$uni}){
                  $self->{unitigs}->{$uni}=undef;
                  delete $self->{unitigs}->{$uni};
                  next;
              }
              #if the unipaths has edges
              my $edges=();
              if(defined $suni->{edges}){
                my $counter=0;
                my $enum=0;
                foreach my $e(@{$suni->{edges}}){
                  $enum++;
                  next if(defined $self->{dunipaths}->{$e->{uid}});
                  push(@{$edges},$e);
                  $counter++;
                }
                #print join(" ","ED",$enum,$counter,$enum-$counter)."\n";
                if($counter > 0){
                      #print Dumper($edges);
                      $suni->{edges}=$edges;
                }else{
                  #the unipaths dont have edges
                  $suni->{edges}=undef;
                }
                #print Dumper($edges);
              }

    }


}


sub print_graph_bcalm{
      my $self= shift;
      my $prefix= shift;
      open(FILE,">".$prefix.".unipath.bcalm.fa") or die "cannot create the $prefix file\n";
      #print Dumper($self->{dunipaths});
      foreach my $uni (keys %{$self->{unitigs}}){
        #we should check how we delete the objects
        next if($uni == "edges");
        my $suni=$self->{unitigs}->{$uni};
        next if(defined $self->{dunipaths}->{$uni});
        #print "U".$suni->get_accession()."\n";
        #print Dumper($suni);
        my @edges=();
        #>3 LN:i:22 KC:i:3 km:f:3.0  L:+:4:- L:+:5:+ L:+:39038:+
        if(defined $suni->{edges}){
          foreach my $e(@{$suni->{edges}}){
            #we skypt the edges that point-out to deleted unipaths;
            next if(defined $self->{dunipaths}->{$e->{uid}});
            #print join(" ","E",$e->{o1},$e->{uid},$e->{o2})."\n";
            my $tmp=join(":","L",$e->{o1},$e->{uid},$e->{o2});
            push(@edges,$tmp);
          }
        }
        print FILE join(" ",">".$suni->get_accession(),"LN:i:".$suni->{LN},"KC:i:".$suni->{KC},"km:f:".$suni->{KM},@edges)."\n";
        print FILE $suni->{sequence}."\n";

      }
      close(FILE);

}





1; #EOM
