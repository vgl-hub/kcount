#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "functions.h" // global functions
#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "threadpool.h"
#include "gfa.h"
#include "sak.h" // swiss army knife
#include "zlib.h"
#include "stream-obj.h"
#include "input-gfa.h"

#include "input.h"
#include "kmer.h"
#include "kcount.h"

void Input::load(UserInputKmap userInput) { // a specialized userInput loading function
    
    this->userInput = userInput;
    
}

void Input::read(short unsigned int mode) { // reads the actual input and performing the tasks
    
    switch (mode) {
            
        case 0: // reads input reads and generates the kmer db
            
        {
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(userInput.kmerLen); // a new empty kmerdb with the specified kmer length
            
            lg.verbose("Loading input sequences");
            unsigned int numFiles = userInput.iReadFileArg.size();
            
            for (unsigned int i = 0; i < numFiles; i++) // for each input read file
                loadKmers(userInput, &kcount, 'r', &i); // specialized function to process reads into kmers as hashes
            
            kcount.hashSegments(); // this is when a fasta is provided
            
            jobWait(threadPool); // ensures that all jobs are done before consolidating the kmerdb
            
            kcount.finalize();
            
            jobWait(threadPool);
            
            lg.verbose("Sequences loaded and hashed");
            
            kcount.hist(); // generates the final histogram
            
            lg.verbose("Histogram generated");
            
            kcount.report(userInput); // output
            
            break;
            
        }
        
        case 1: // reads an existing kmerdb
            
        {
            
            std::ifstream file;
            
            file.open(userInput.iSeqFileArg + "/.index"); // reads the kmer length from the index file
            std::string line;
            
            getline(file, line);
            file.close();
            
            short unsigned int k = stoi(line);
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(k); // a new empty kmerdb with the specified kmer length
            
            lg.verbose("Kmer object generated");
            
            kcount.load(userInput); // loads kmers into the new kmerdb
            
            kcount.report(userInput); // output
            
            break;
            
        }
            
        case 2: // union of multiple kmerdbs
            
        {
            
            std::ifstream file;
            
            lg.verbose("Merging input databases");
            unsigned int numFiles = userInput.iReadFileArg.size(); // number of input kmerdbs
            
            short unsigned int k = 0;
            
            for (unsigned int i = 0; i < numFiles; i++) {  // reads the kmer length from the index files checking consistency between kmerdbs
                
                file.open(userInput.iReadFileArg[i] + "/.index");
                std::string line;
                
                getline(file, line);
                file.close();
                
                if (k == 0)
                    k = stoi(line);
                
                if (k != stoi(line)) {
                    fprintf(stderr, "Cannot merge databases with different kmer length\n");
                    exit(1);
                }
                
            }
            
            if (k == 0 || k > 32) {
                fprintf(stderr, "Invalid kmer length\n");
                exit(1);
            }
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(k); // a new empty kmerdb with the specified kmer length
            
            lg.verbose("Kmer object generated. Merging.");
            
            kcount.kunion(userInput); // union set
            
            kcount.report(userInput); // output
            
            break;
            
        }
            
        default:
            
            fprintf(stderr, "Invalid mode\n");
            exit(1);
        
    }
    
}
