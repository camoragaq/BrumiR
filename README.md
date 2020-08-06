# BrumiR

A de bruijn graph-based tool for miRNA discovery.

## Running Brumir

```

perl brumir.pl -a test/sRNA-seq.human.trim.fa.gz -p prefix

```

### BrumiR outputs

| File  |  Description  |   
|-------|---------------|
| prefix.brumir.candidate_miRNA.fasta   |  fasta file with all the candidates with their KM and KC values respectively. |
|  prefix.brumir.other_sequences.txt |  fasta file with all long sequences expressed in the sample, their are putative long non-coding RNA. |
| prefix.brumir.RFAM_HITS.txt | table with a list of putative tRNAs or rRNAs presents in RFAM database. |

 

## Running BrumiR2reference

```bash
perl brumir2reference.pl -a prefix.candidate_miRNA.fasta -b test/chr1-20M-50M.human.fna -p prefix2ref
```
the file *.passfilter.txt contains the miRNAs that have a valid precursor sequence.

## Convert BrumiR output to GFA (Bandage)
```python aux_scripts/convertToGFA.py h.test.unipath.bcalm.fa h.test.unipath.bcalm.gfa 18```


# Getting the latest source code
## Instructions
It is recommended to use/download the latest binary release (Linux or Mac) from : https://github.com/camoragaq/BrumiR/releases

## Building BrumiR from souce

To compile BrumiR run the following command:
```bash
#fetch BrumiR 
git clone https://github.com/camoragaq/BrumiR.git brumir
cd brumir
make all
```

### Dependencies

#### BCALM (v2.2.2)
The BCALM binary can be donloaded from: 
##### Linux binaries
```
wget https://github.com/GATB/bcalm/releases/download/v2.2.2/bcalm-binaries-v2.2.2-Linux.tar.gz
```
##### Mac binaries
```
wget https://github.com/GATB/bcalm/releases/download/v2.2.2/bcalm-binaries-v2.2.2-Mac.tar.gz
```
The bcalm binary should be copied to BrumiR-dir/bin


#### RNAfold (v2.2.14)
download the ViennaRNA package and follow the instructions to compile the code.
```bash
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz
tar -zxvf ViennaRNA-2.4.14.tar.gz
./configure
make
make install
```

The RNAfold binary should be copied to BrumiR-dir/bin
