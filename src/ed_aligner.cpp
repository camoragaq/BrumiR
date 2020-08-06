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



#include "edlib.h"

using namespace std;

//struct for seqs
typedef  struct seq{
    int id=-1;
    string name="";
    string seq="";//contigs are
} seq;


int read_fasta(string filename, vector<seq> &seqs);


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
    int option;
    int kArg = -1;

    // If "STD" or "EXT", cigar string will be printed. if "NICE" nice representation
    // of alignment will be printed.
    char alignmentFormat[16] = "NICE";

    bool invalidOption = false;
    while ((option = getopt(argc, argv, "m:n:k:f:spl")) >= 0) {
        switch (option) {
        case 'm': strcpy(mode, optarg); break;
        case 'n': numBestSeqs = atoi(optarg); break;
        case 'k': kArg = atoi(optarg); break;
        case 'f': strcpy(alignmentFormat, optarg); break;
        case 's': silent = true; break;
        case 'p': findAlignment = true; break;
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
        fprintf(stderr, "\t-p  If specified, alignment path will be found and printed. "
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
        queriesTotalLength += queries[i].seq.length();
	//printf("id=%d name=%s seq1=%s\n",i, queries[i].name.c_str() ,queries[i].seq.c_str());
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
    for(auto t=0; t < targets.size(); t++){		
    const char* target = targets[t].seq.c_str();
    const char* ntarget = targets[t].name.c_str();
    int targetLength = targets[t].seq.length();
    printf("Read target name=%s id=%d,L=%d residues.\n",ntarget,t, targetLength);


    // ----------------------------- MAIN CALCULATION ----------------------------- //
    int* scores = new int[numQueries];
    int** endLocations = new int*[numQueries];
    int** startLocations = new int*[numQueries];
    int* numLocations = new int[numQueries];
    priority_queue<int> bestScores; // Contains numBestSeqs best scores
    int k = kArg;
    unsigned char* alignment = NULL; int alignmentLength;

   /* if (!findAlignment || silent) {
        printf("0/%d", numQueries);
        fflush(stdout);
    }
   */
    for (int i = 0; i < numQueries; i++) {

        const char* query = queries[i].seq.c_str();
        int queryLength = queries[i].seq.length();
        // Calculate score
        EdlibAlignResult result = edlibAlign(query, queryLength, target, targetLength,
                                             edlibNewAlignConfig(k, modeCode, alignTask, NULL, 0));
        scores[i] = result.editDistance;
        endLocations[i] = result.endLocations;
        startLocations[i] = result.startLocations;
        numLocations[i] = result.numLocations;
        alignment = result.alignment;
        alignmentLength = result.alignmentLength;

        // If we want only numBestSeqs best sequences, update best scores 
        // and adjust k to largest score.
        if (numBestSeqs > 0) {
            if (scores[i] >= 0) {
                bestScores.push(scores[i]);
                if ((int) bestScores.size() > numBestSeqs) {
                    bestScores.pop();
                }
                if ((int) bestScores.size() == numBestSeqs) {
                    k = bestScores.top() - 1;
                    if (kArg >= 0 && kArg < k)
                        k = kArg;
                }
            }
        }
        
        if (!findAlignment || silent) {
           //printf("\r%d/%d", i + 1, numQueries);
            //fflush(stdout);
        } else {
            // Print alignment if it was found, use first position
            if (alignment) {
                printf("\n");
                //printf("Query #%d (%d residues): score = %d\n", i, queryLength, scores[i]);
                printf("#%s %s: %d  %d\n", ntarget, queries[i].name.c_str(), scores[i], numLocations[i]);
                if (!strcmp(alignmentFormat, "NICE")) {
                    printAlignment(query, target, alignment, alignmentLength,
                                   *(endLocations[i]), modeCode);
                } else {
                    printf("Cigar:\n");
                    EdlibCigarFormat cigarFormat = !strcmp(alignmentFormat, "CIG_STD") ?
                        EDLIB_CIGAR_STANDARD : EDLIB_CIGAR_EXTENDED;
                    char* cigar =edlibAlignmentToCigar(alignment, alignmentLength, cigarFormat);
                    if (cigar) {
                        printf("%s\n", cigar);
                        free(cigar);
                    } else {
                        printf("Error while printing cigar!\n");
                    }
                }
            }
        }

        if (alignment)
            free(alignment);
    }

    if (!silent && !findAlignment) {
        int scoreLimit = -1; // Only scores <= then scoreLimit will be printed (we consider -1 as infinity)
        printf("\n");

        if (bestScores.size() > 0) {
            printf("%d best scores:\n", (int)bestScores.size());
            scoreLimit = bestScores.top();
        } else {
            printf("Scores:\n");
        }

        for (int i = 0; i < numQueries; i++) {
            if (scores[i] > -1 && (scoreLimit == -1 || scores[i] <= scoreLimit)) {
                //printf("#%d: %d  %d", i, scores[i], numLocations[i]);
                printf("#%s %s: %d  %d", ntarget, queries[i].name.c_str(), scores[i], numLocations[i]);
                if (numLocations[i] > 0) {
                    printf("  [");
                    for (int j = 0; j < numLocations[i]; j++) {
                        printf(" (");
                        if (startLocations[i]) {
                            printf("%d", *(startLocations[i] + j));
                        } else {
                            printf("0");
                        }
                        printf(", %d)", *(endLocations[i] + j));
                    }
                    printf(" ]");
                }
                printf("\n");
            }
        }

    }

   // ---------------------------------------------------------------------------- //

    // Free allocated space
    for (int i = 0; i < numQueries; i++) {
        free(endLocations[i]);
        if (startLocations[i]) free(startLocations[i]);
    }
    delete[] endLocations;
    delete[] startLocations;
    delete[] numLocations;
    delete[] scores;
    } //closing the free 

    clock_t finish = clock();
    double cpuTime = ((double)(finish-start))/CLOCKS_PER_SEC;
    printf("\nCpu time of searching: %lf\n", cpuTime);
   
    return 0;
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


