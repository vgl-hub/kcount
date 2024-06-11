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

#include <parallel-hashmap/phmap.h>

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
            
            Kmap<UserInputKmap, uint8_t, uint32_t> kcount(userInput); // a new empty kmerdb with the specified kmer length
            
            lg.verbose("Loading input sequences");
            unsigned int numFiles = userInput.inReads.size();
            
            for (unsigned int i = 0; i < numFiles; ++i) // for each input read file
                loadKmers(userInput, &kcount, 'r', i); // specialized function to process reads into kmers as hashes
            
            lg.verbose("Reads loaded.");
            kcount.finalize();
            
            kcount.loadHighCopyKmers();
            kcount.report(userInput); // output
            kcount.cleanup(); // delete tmp files
            break;
            
        }
        
        case 1: // reads an existing kmerdb
            
        {
            
            std::ifstream file;
            
            file.open(userInput.inSequence + "/.index"); // reads the kmer length from the index file
            std::string line;
            
            getline(file, line);
            file.close();
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(userInput); // a new empty kmerdb with the specified kmer length
            
            lg.verbose("Kmer object generated");
            
            kcount.load(userInput); // loads kmers into the new kmerdb
            
            kcount.stats();
            
            kcount.report(userInput); // output
            
            break;
            
        }
            
        case 2: // union of multiple kmerdbs
            
        {
            
            std::ifstream file;
            
            lg.verbose("Merging input databases");
            unsigned int numFiles = userInput.inReads.size(); // number of input kmerdbs
            
            short unsigned int k = 0;
            
            for (unsigned int i = 0; i < numFiles; i++) {  // reads the kmer length from the index files checking consistency between kmerdbs
                
                file.open(userInput.inReads[i] + "/.index");
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
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(userInput); // a new empty kmerdb with the specified kmer length
            
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
