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

void Input::load(UserInputKmap userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(short unsigned int mode) {
    
    switch (mode) {
            
        case 0:
            
        {
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(userInput.kmerLen);
            
            lg.verbose("Loading input sequences");
            unsigned int numFiles = userInput.iReadFileArg.size();
            
            for (unsigned int i = 0; i < numFiles; i++)
                loadKmers(userInput, &kcount, 'r', &i);
            
            kcount.hashSegments();
            
            jobWait(threadPool);
            
            kcount.finalize();
            
            jobWait(threadPool);
            
            lg.verbose("Sequences loaded and hashed");
            
            kcount.hist();
            
            lg.verbose("Histogram generated");
            
            kcount.report(userInput);
            
            break;
            
        }
        
        case 1:
            
        {
            
            std::ifstream file;
            
            file.open(userInput.iSeqFileArg + "/.index");
            std::string line;
            
            getline(file, line);
            file.close();
            
            short unsigned int k = stoi(line);
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(k);
            
            lg.verbose("Kmer object generated");
            
            kcount.load(userInput);
            
            kcount.report(userInput);
            
            break;
            
        }
            
        case 2:
            
        {
            
            std::ifstream file;
            
            lg.verbose("Merging input databases");
            unsigned int numFiles = userInput.iReadFileArg.size();
            
            short unsigned int k = 0;
            
            for (unsigned int i = 0; i < numFiles; i++) {
                
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
            
            Kmap<UserInputKmap, uint32_t, uint64_t> kcount(k);
            
            lg.verbose("Kmer object generated. Merging.");
            
            kcount.kunion(userInput);
            
            kcount.report(userInput);
            
            break;
            
        }
            
        default:
            
            fprintf(stderr, "Invalid mode\n");
            exit(1);
        
    }
    
}
