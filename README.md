<p align="center"><img src="img/logo.png"
alt="BrumiR" width="400" height="264" border="0" /></p>

# BrumiR

A de Bruijn graph-based tool for miRNA discovery.

Description
============

miRNAs are small non-coding RNAs that have become key players at the regulation level. In the last decade, with the increase and accessibility of high-throughput sequencing technologies, different methods have been developed to identify these molecules, but most of them are genome-based methods or are methods based on conservation information. However, when high quality reference genomes are not available, our possibilities are considerably reduced. In this context, we developed BrumiR, a de novo algorithm based on a de Bruijn approach, which is able to identify miRNAs directly and exclusively from sRNA-seq data. We tested BrumiR on different datasets (simulated and real sRNA-seq) from the animal and plant kingdoms, and we compared the results with the state-of-the-art tools in this field. The results of BrumiR exceeded or were comparable to those of the existing methods. Moreover, BrumiR is an ultra-fast algorithm, 20X faster than the state-of-the-art tools, enabling the analysis of a large number of experiments. Futhermore, BrumiR is very easy to use and provides additional tools to explore the results. It also identifies other small RNAs in order to maximize the biological insight. BrumiR presents a new and versatile method that implements novel algorithmic ideas for the study of miRNAs that complements and extends the currently existing approaches. 

## Running Brumir

```

perl brumir.pl -a test/sRNA-seq.human.trim.fa.gz -p prefix

```

### BrumiR outputs

| File  |  Description  |   
|-------|---------------|
| prefix.brumir.candidate_miRNA.fasta   |  fasta file with all the candidates with their KM and KC values respectively. |
|  prefix.brumir.other_sequences.txt |  asta file with all long sequences expressed in the sample, they are putative long non-coding RNAs. |
| prefix.brumir.RFAM_HITS.txt | table with a list of putative tRNAs or rRNAs present in the RFAM database. |

 

## Running BrumiR2reference

```bash
perl brumir2reference.pl -a prefix.candidate_miRNA.fasta -b test/chr1-20M-50M.human.fna -p prefix2ref
```
the file *.passfilter.txt contains the miRNAs that have a valid precursor sequence.

## Convert BrumiR output to GFA (Bandage)
```python aux_scripts/convertToGFA.py h.test.unipath.bcalm.fa h.test.unipath.bcalm.gfa 18```

## Getting the latest source code
## Instructions
It is recommended to use/download the latest binary release (Linux or Mac) from : https://github.com/camoragaq/BrumiR/releases


### Containers
To facilitate the execution of BrumiR, we provide docker/singularity containers.
BrumiR images are hosted on [Dockerhub](https://hub.docker.com/repository/docker/camoragaq/brumir) and can be downloaded with the command:

```
docker pull camoragaq/brumir:v2.0
```

Alternatively, using singularity:

```
export TMPDIR=/tmp
singularity pull docker://camoragaq/brumir:v2.0
```

#### Run BrumiR using singularity
```
#using singularity
CONTAINER=/camoragaq/brumir_v2.0.sif


#run BrumiR with singularity exec
singularity exec brumir_v2.0.sif perl brumir.pl --help



## Building BrumiR from source

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


