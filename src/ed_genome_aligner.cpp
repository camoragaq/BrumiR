
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <climits>
#include <queue>
#include <map>
#include <unordered_map>
#include <iterator>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>


#include "edlib.h"

using namespace std;

//struct for seqs
typedef  struct seq{
    int id=-1;
    string name="";
    string seq="";//contigs are
    string rseq="";
} seq;

//save the match of miRNAs in the genome sequence
typedef struct mihit{
  int qid;// internal microRNA id
  string qname;//name of the miRNA
  int rid;//internal reference id
  int qstart;//miStart
  int qend;//Mi end
  int tstart;//ref start
  int tend;//ref name
  int edit_distance;//actual edit distance
  bool strand;//0 means forward, 1 means reverse
  int alnlen;//length of alignment
  int qlen;//miRNA lenght
  int cov;//aln coverage of the miRNA
}mihit;

string revcomp(string kmer);

void map_in_genome_window(string wgenome, int wstart, vector<seq> &queries, vector<mihit> &mihits, int k, bool log,int t);

int read_fasta(string filename, vector<seq> &seqs);

void print_hit(mihit &mh);

//function that get the sequence of the potencial precursors
void get_potential_precursors(vector<mihit> &mihits, vector<seq> &qry, vector<seq> &ref, int bestn);
//function that get the sequence of the potencial precursors for plants
void get_potential_precursors_plants(vector<mihit> &mihits, vector<seq> &qry, vector<seq> &ref, int bestn);

void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode);

// For debugging
void printSeq(const vector<char> &seq) {
    for (int i = 0; i < (int) seq.size(); i++)
        printf("%d ", seq[i]);
    printf("\n");
}

int main(int argc, char * const argv[]) {

    //----------------------------- PARSE COMMAND LINE ------------------------//
    // If true, there will be no output.
    bool silent = false;
    // Alignment mode.
    char mode[16] = "NW";
    // How many best sequences (those with smallest score) do we want.
    // If 0, then we want them all.
    int numBestSeqs = 0;
    bool findAlignment = false;
    bool findStartLocations = false;
    bool plantmode = false;
    int option;
    int kArg = -1;

    // If "STD" or "EXT", cigar string will be printed. if "NICE" nice representation
    // of alignment will be printed.
    char alignmentFormat[16] = "NICE";

    bool invalidOption = false;
    while ((option = getopt(argc, argv, "m:n:k:f:spla")) >= 0) {
        switch (option) {
        case 'm': strcpy(mode, optarg); break;
        case 'n': numBestSeqs = atoi(optarg); break;
        case 'k': kArg = atoi(optarg); break;
        case 'f': strcpy(alignmentFormat, optarg); break;
        case 's': silent = true; break;
        case 'p': plantmode = true; break;
        case 'a': findAlignment = true; break;
        case 'l': findStartLocations = true; break;
        default: invalidOption = true;
        }
    }
    if (optind + 2 != argc || invalidOption) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: %s [options...] <queries.fasta> <target.fasta>\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "\t-s  If specified, there will be no score or alignment output (silent mode).\n");
        fprintf(stderr, "\t-m HW|NW|SHW  Alignment mode that will be used. [default: NW]\n");
        fprintf(stderr, "\t-n N  Score will be calculated only for N best sequences (best = with smallest score)."
                " If N = 0 then all sequences will be calculated."
                " Specifying small N can make total calculation much faster. [default: 0]\n");
        fprintf(stderr, "\t-k K  Sequences with score > K will be discarded."
                " Smaller k, faster calculation. If -1, no sequences will be discarded. [default: -1]\n");
        fprintf(stderr, "\t-p  If specified, active the plante mode of the genome based mapping.\n ");
        fprintf(stderr, "\t-a  If specified, alignment path will be found and printed. "
                        "This may significantly slow down the calculation.\n");
        fprintf(stderr, "\t-l  If specified, start locations will be found and printed. "
                "Each start location corresponds to one end location. This may somewhat slow down "
                "the calculation, but is still faster then finding alignment path and does not consume "
                "any extra memory.\n");
        fprintf(stderr, "\t-f NICE|CIG_STD|CIG_EXT  Format that will be used to print alignment path,"
                " can be used only with -p. NICE will give visually attractive format, CIG_STD will "
                " give standard cigar format and CIG_EXT will give extended cigar format. [default: NICE]\n");
        return 1;
    }
    //-------------------------------------------------------------------------//

    if (strcmp(alignmentFormat, "NICE") && strcmp(alignmentFormat, "CIG_STD") &&
        strcmp(alignmentFormat, "CIG_EXT")) {
        printf("Invalid alignment path format (-f)!\n");
        return 1;
    }

    EdlibAlignMode modeCode;
    if (!strcmp(mode, "SHW"))
        modeCode = EDLIB_MODE_SHW;
    else if (!strcmp(mode, "HW"))
        modeCode = EDLIB_MODE_HW;
    else if (!strcmp(mode, "NW"))
        modeCode = EDLIB_MODE_NW;
    else {
        printf("Invalid mode (-m)!\n");
        return 1;
    }
    printf("Using %s alignment mode.\n", mode);

    EdlibAlignTask alignTask = EDLIB_TASK_DISTANCE;
    if (findStartLocations) alignTask = EDLIB_TASK_LOC;
    if (findAlignment) alignTask = EDLIB_TASK_PATH;


    int readResult;
    // Read queries
    char* queriesFilepath = argv[optind];
    vector <seq> queries;

    printf("Reading queries...\n");
    readResult = read_fasta(queriesFilepath, queries);
    if (readResult) {
        printf("Error: There is no file with name %s\n", queriesFilepath);
        return 1;
    }



    int numQueries = queries.size();
    printf("NQ=%d\n",numQueries);
    int queriesTotalLength = 0;
    for (int i = 0; i < numQueries; i++) {
        queries[i].rseq=revcomp(queries[i].seq);
        queriesTotalLength += queries[i].seq.length();
	       printf("id=%d name=%s seq1=%s rseq1=%s\n",i, queries[i].name.c_str() ,queries[i].seq.c_str(),queries[i].rseq.c_str());
    }
    printf("Read %d queries, %d residues total.\n", numQueries, queriesTotalLength);

    // Read target
    char* targetFilepath = argv[optind+1];
    printf("Reading target fasta file...\n");
   vector<seq> targets;
    readResult = read_fasta(targetFilepath, targets);
    if (readResult) {
        printf("Error: There is no file with name %s\n", targetFilepath);
        return 1;
    }

    clock_t start = clock();
    printf("\nComparing queries to targets...\n");
    //for quering all agains all
    vector<mihit> mihits;
    int window=200;
    for(auto t=0; t < targets.size(); t++){
      //we split target in windows of length 500bp
      auto starget=targets[t].seq;
        printf("Scanning seq:%s %d\n",targets[t].name.c_str(),targets[t].seq.length());
      for(auto j=0; j<starget.length()-window; j+=window) {
          auto wgenome = starget.substr(j, window);
          //we map the microRNA in the window
          map_in_genome_window(wgenome, j, queries, mihits, 2, false, t);
      }
  /*  const char* target = targets[t].seq.c_str();
    const char* ntarget = targets[t].name.c_str();
    int targetLength = targets[t].seq.length();*/
    //printf("Read target name=%s id=%d,L=%d residues.\n",ntarget,t, targetLength);

    } //closing the free
    printf("Total hits stored %d\n",mihits.size());
    //for(auto h:mihits)
    //sort(mihits.begin(),mihits.end())
  //  for(auto mh:mihits)
  //      print_hit(mh);
    //Step 2, we sort the array of positions by minimizer,seq and pos
    sort( mihits.begin( ), mihits.end( ), [ ]( const mihit a, const mihit b)
    {
        //return tie(a.qid,a.edit_distance,a.rid,a.tstart) < tie(b.qid,b.edit_distance,b.rid,b.tstart);
        //return tie(a.qid,a.edit_distance,a.cov*-1) < tie(a.qid,b.edit_distance,b.cov*-1);
        return (a.qid < b.qid) ||
       ((a.qid == b.qid) && (a.edit_distance < b.edit_distance)) ||
       ((a.qid == b.qid) && (a.edit_distance == b.edit_distance) &&
          (a.cov > b.cov));


    });

    //we select the best N hits per miRNA location
    int MAX_HITS_PER_MIRNA=100;
    printf("Selecting the best %d per miRNA\n",MAX_HITS_PER_MIRNA);
    vector<mihit> best_mihits;
    unordered_map<int, int> hitcounter;
    //we innit the hitcounter hash
    for( auto mi:queries){
          hitcounter[mi.id]=0;
    }
    //we select the best MAX_HITS_PER_MIRNA
    for(auto mh:mihits){
        if(hitcounter[mh.qid] < MAX_HITS_PER_MIRNA){
                    best_mihits.push_back(mh);
                    hitcounter[mh.qid]++;
          }else{
                    hitcounter[mh.qid]++;
          }
    }
    //we print the total HITS per miRNA
    for(auto cc:hitcounter)
          printf("H_TOTAL_COUNT: id=%d total=%d\n",cc.first,cc.second);

    //we print the best hit selected
    for(auto mh:best_mihits)
        print_hit(mh);

    //get the sequece from the hits
    if(plantmode){
        //plante mode
        //get_potential_precursors_plants(mihits,queries,targets,10);
        get_potential_precursors_plants(best_mihits,queries,targets,10);
    }else{
          //default mode for vertebrates genomes
          //get_potential_precursors(mihits,queries,targets,10);
          get_potential_precursors(best_mihits,queries,targets,10);
    }
    clock_t finish = clock();
    double cpuTime = ((double)(finish-start))/CLOCKS_PER_SEC;
    printf("\nCpu time of searching: %lf\n", cpuTime);

    return 0;
}

//global array for  to compute revcomp
char const comp_tab2[] = {
        0,   1,       2,       3,       4,   5,       6,       7,       8,   9,  10,  11,      12,  13,  14,  15,
        16,  17,  18,  19,      20,  21,  22,  23,      24,  25,  26,  27,      28,  29,  30,  31,
        32,  33,  34,  35,      36,  37,  38,  39,      40,  41,  42,  43,      44,  45,  46,  47,
        48,  49,  50,  51,      52,  53,  54,  55,      56,  57,  58,  59,      60,  61,  62,  63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,      92,  93,  94,  95,
        64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

//heng li revcomp this is for printing only

string revcomp(string kmer) {
    int c0, c1;
    int klen=kmer.length();
    for (int i = 0; i < klen>>1; ++i) { // reverse complement sequence
        c0 = comp_tab2[(int)kmer[i]];
        c1 = comp_tab2[(int)kmer[klen - 1 - i]];
        kmer[i] = c1;
        kmer[klen - 1 - i] = c0;
    }
    if (klen & 1) // complement the remaining base
        kmer[klen>>1] = comp_tab2[(int)kmer[klen>>1]];
    return kmer;
}

void print_hit(mihit &mh){
    printf("HIT Qid=%d Rid=%d rs=%d re=%d ed=%d aln=%d  qlen=%d cov=%d strand=%d\n",mh.qid,mh.rid,mh.tstart,mh.tend,mh.edit_distance,mh.alnlen,mh.qlen,mh.cov,mh.strand);
}


//mirDeep2 take the sequences like this
//todo: check RNAfold command of mirDeep2
//ViennaRNA-2.4.9/RNA/bin/RNAfold -j3  --noPS -i ara.hits.fa --outfile=ara.hits.rnafold, three is the number of cores
//command from mirDeep2 RNAfold < $dir_tmp/precursors.fa -noPS > $dir_tmp/precursors.str 2>>error_${time}.log
//at this point the array of hits is sorted by miRNA and we want to create the fasta file for the precursor
void get_potential_precursors(vector<mihit> &mihits, vector<seq> &qry, vector<seq> &ref, int bestn){
    for (auto mh:mihits) {
               //print_hit(mh);
               //we print the foward sequence
                int excise_beg= mh.tstart-70;
                int excise_end= mh.tend+20;
                if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                //printf("hsp_%d_%d_%d_%d %s\n",mh.qid,mh.rid,excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());
                printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());

                //we print the rev sequence
                excise_beg= mh.tstart-20;
                excise_end= mh.tend+70;
                if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                //printf("hsp_%d_%d_%d_%d %s\n",mh.qid,mh.rid,excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());
                printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());

    }
}


//print precursor sequence for plants
void get_potential_precursors_plants(vector<mihit> &mihits, vector<seq> &qry, vector<seq> &ref, int bestn){
    for (auto mh:mihits) {
               //print_hit(mh);
               //create a precursor sequence of ~110bp
               //we print the foward sequence
                int excise_beg= mh.tstart-70;
                int excise_end= mh.tend+20;
                if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());
                //we print the rev sequence
                excise_beg= mh.tstart-20;
                excise_end= mh.tend+70;
                if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());

                //create a precursor sequence of ~150bp
                  excise_beg= mh.tstart-120;
                  excise_end= mh.tend+20;
                 if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                 printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());
                 //we print the rev sequence
                 excise_beg= mh.tstart-20;
                 excise_end= mh.tend+120;
                 if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                 printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());

                 //create a precursor sequence of ~200bp
                   excise_beg= mh.tstart-160;
                   excise_end= mh.tend+20;
                  if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                  printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());
                  //we print the rev sequence
                  excise_beg= mh.tstart-20;
                  excise_end= mh.tend+160;
                  if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                  printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());


                  //create a precursor sequence of ~250bp
                    excise_beg= mh.tstart-210;
                    excise_end= mh.tend+20;
                   if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                   printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());
                   //we print the rev sequence
                   excise_beg= mh.tstart-20;
                   excise_end= mh.tend+210;
                   if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                   printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());

                   //create a precursor sequence of ~300bp
                     excise_beg= mh.tstart-260;
                     excise_end= mh.tend+20;
                    if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                    printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());
                    //we print the rev sequence
                    excise_beg= mh.tstart-20;
                    excise_end= mh.tend+260;
                    if(excise_beg >=0 and excise_end <=ref[mh.rid].seq.length())
                    printf("hsp::%s::%s::%d::%d %s\n",mh.qname.c_str(),ref[mh.rid].name.c_str(),excise_beg,excise_end,ref[mh.rid].seq.substr(excise_beg,abs(excise_beg-excise_end)).c_str());

    }
}



/*
 * https://github.com/rajewsky-lab/mirdeep2/blob/master/src/excise_precursors.pl
 * lines: 251-261
 * 	  #else excise to sequences, corresponding to the read stack being the mature sequence
		#in the 5' arm or the 3' arm of the precursor hairpin
		my $excise_beg=$db_beg-70;
		my $excise_end=$db_end+20;

		#print out in fasta format
		print_positions($db,\$strand,$db_seq,\$db_lng,\$excise_beg,\$excise_end);
		$excise_beg=$db_beg-20;
		$excise_end=$db_end+70;
		print_positions($db,\$strand,$db_seq,\$db_lng,\$excise_beg,\$excise_end);
		#the most 3' position that has yet been excised
		$db_limit=$excise_end;
 */





void map_in_genome_window(string wgenome, int wstart, vector<seq> &queries, vector<mihit> &mihits,int k, bool log, int t){

  EdlibAlignTask alignTask = EDLIB_TASK_LOC;
  //if (findStartLocations) alignTask = EDLIB_TASK_LOC;
  if (log) alignTask = EDLIB_TASK_PATH;
    EdlibAlignMode modeCode;
    modeCode=EDLIB_MODE_HW;
  const char* target  = wgenome.c_str();
  //hay que hacer lo mismo para el reverse
  int kmersize=12;//seed size miRNAs
  mihit mh;
  //indexando los 10-mers de la ventana
  unordered_map<string, bool> hashqq;
  for(int j=0; j< wgenome.length()-kmersize; j++)
        hashqq[wgenome.substr(j,kmersize)]=true;
//printf("%d %s\n",queries.size(),target);
  for (int i = 0; i < queries.size(); ++i) {
        //printf("%s %s\n",queries[i].seq.c_str(),target);
      //char* query = queries[i].seq.c_str();
      //int queryLength = queries[i].seq.length();
      bool aln=false;
      for(int k=0; k<queries[i].seq.length()-kmersize; k++)
          if(hashqq.find(queries[i].seq.substr(k,kmersize)) !=hashqq.end() ){
            aln=true;
            break;
          }

      if(aln){
      // Calculate score
      EdlibAlignResult result = edlibAlign(queries[i].seq.c_str(), queries[i].seq.length(), target, wgenome.length(),
                                           edlibNewAlignConfig(-1, modeCode, alignTask, NULL, 0));

      if (result.status == EDLIB_STATUS_OK){
          if(result.editDistance <= k) {
              /*printf("HIT %s %s %d %d %d %d %d %d\n", queries[i].name.c_str(), rname.c_str(),
                     wstart, result.editDistance, result.alignmentLength,
                     result.startLocations[0], result.endLocations[0], 0);*/
              //mh.qname = queries[i].name;
              mh.qid=i;
              mh.qname=queries[i].name;
              mh.rid=t;
              //mh.rname = rname;
              mh.tstart = wstart + result.startLocations[0];//pos for the target
              mh.tend = wstart + result.endLocations[0];
              mh.qstart = result.startLocations[0];
              mh.qend = result.endLocations[0];
              mh.edit_distance = result.editDistance;
              mh.strand = 0;//strand foward for the miRNA
              //mh.alnlen = result.alignmentLength;
              mh.alnlen = abs(result.startLocations[0]-result.endLocations[0])+1;//is 0-based
              mh.qlen=queries[i].seq.length();
              mh.cov=((float)mh.alnlen/mh.qlen)*100;
              //print_hit(mh);
              mihits.push_back(mh);//we store the hit in the vector
          }
        }
      edlibFreeAlignResult(result);
      }
      //reverse strand
      aln=false;
      for(int k=0; k<queries[i].rseq.length()-kmersize; k++)
          if(hashqq.find(queries[i].rseq.substr(k,kmersize)) !=hashqq.end() ){
            aln=true;
            break;
          }

      if(aln){
      // Calculate score
      EdlibAlignResult result = edlibAlign(queries[i].rseq.c_str(), queries[i].rseq.length(), target, wgenome.length(),
                                           edlibNewAlignConfig(-1, modeCode, alignTask, NULL, 0));

      if (result.status == EDLIB_STATUS_OK){
          if(result.editDistance <= k) {
/*
              printf("HIT %s %s %d %d %d %d %d %d\n", queries[i].name.c_str(), rname.c_str(),
                     wstart, result.editDistance, result.alignmentLength,
                     result.startLocations[0], result.endLocations[0], 1);
*/
              //mh.qname = queries[i].name;
              mh.qid=i;
              mh.qname=queries[i].name;
              mh.rid=t;
              //mh.rname = rname;
              mh.tstart = wstart + result.startLocations[0];//pos for the target
              mh.tend = wstart + result.endLocations[0];
              mh.qstart = result.startLocations[0];
              mh.qend = result.endLocations[0];
              mh.edit_distance = result.editDistance;
              mh.strand = 1;//strand reverse of the miRNA
              mh.alnlen = abs(result.startLocations[0]-result.endLocations[0])+1;
              mh.qlen=queries[i].seq.length();
              mh.cov=((float)mh.alnlen/mh.qlen)*100;
              //print_hit(mh);
              mihits.push_back(mh);//we store the hit in the vector
          }
      }
      edlibFreeAlignResult(result);
      }

  }

}



void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        printf("T: ");
        int startTIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                printf("-");
            else
                printf("%c", target[++tIdx]);
            if (j == start)
                startTIdx = tIdx;
        }
        printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);

        // match / mismatch
        printf("   ");
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        printf("\n");

        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                printf("-");
            else
                printf("%c", query[++qIdx]);
            if (j == start)
                startQIdx = qIdx;
        }
        printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
    }
}

//Construct contigs from FASTA file
int read_fasta(string filename, vector<seq> &seqs){
    std::ifstream input(filename);
    if(!input.good()){
        std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
        return 1;
    }
    int id=0;
    std::string line, name, content;
    while(getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << ":" << content << std::endl;
                //cout << name << ":" <<id<<":"<< content << std::endl;
                //todo: improve the parsing of contig name for the moment we expect and space
                seq tmp;
                tmp.id=id;
                tmp.name=name.find(' ') ? name.substr(0, name.find(' ')):name.substr(0, name.find('\n'));
                tmp.seq=content;
		seqs.push_back(tmp);

                id++;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);//pick name before the end;
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }

    if( !name.empty() ){ // Print out what we read from the last entry
	 seq tmp;
         tmp.id=id;
         tmp.name=name.substr(0, name.find(' '));
         tmp.seq=content;
	 seqs.push_back(tmp);
        //cout << name << ":" <<id<<":"<< content << std::endl;
    }
	return 0;
}
