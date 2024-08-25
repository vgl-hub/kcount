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

void Input::load(UserInputKcount userInput) { // a specialized userInput loading function
    
    this->userInput = userInput;
    
}

void Input::loadDB() {
    
    if (userInput.kmerDB.size() == 1){
        userInput.prefix = userInput.kmerDB[0]; // access database
        std::ifstream file;
        file.open(userInput.prefix + "/.index"); // update kmer length
        std::string line;
        getline(file, line);
        file.close();
        userInput.kmerLen = stoi(line);
        kLen = userInput.kmerLen;
//        if (kLen < 31)
//            kPrefixLen = kLen;
    }else if (userInput.kmerDB.size() > 1) {
        fprintf(stderr, "More than one DB provided. Merge them first. Exiting.\n");
        exit(EXIT_FAILURE);
    }else{
        fprintf(stderr, "Cannot load DB input. Exiting.\n");
        exit(EXIT_FAILURE);
    }
    
}

void Input::read() { // reads the actual input and performing the tasks
    
    if (userInput.outFile.find(".kc") != std::string::npos)
        userInput.prefix = userInput.outFile;
    
    if (userInput.prefix != ".")
        make_dir(userInput.prefix.c_str());
    
    switch (userInput.mode) {
            
        case 0: { // reads input reads and generates the kmer db
            
            KDB kcount(userInput); // a new empty kmerdb with the specified kmer length
            
            if (userInput.inReads.size() > 0) {
                
                lg.verbose("Loading input sequences");
                unsigned int numFiles = userInput.inReads.size();
                
                for (unsigned int i = 0; i < numFiles; ++i) // for each input read file
                    loadKmers(userInput, kcount, 'r', i); // specialized function to process reads into kmers as hashes
                
                lg.verbose("Reads loaded.");
                kcount.finalize();
            }else{
                fprintf(stderr, "Reads not provided. Exiting.\n");
                exit(EXIT_FAILURE);
            }
            kcount.report(); // output
            kcount.cleanup(); // delete tmp files
            break;
        }
        
        case 1: { // reads an existing kmerdb
            
            std::ifstream file;
            file.open(userInput.inSequence + "/.index"); // reads the kmer length from the index file
            std::string line;
            getline(file, line);
            file.close();
            
            KDB kcount(userInput); // a new empty kmerdb with the specified kmer length
            lg.verbose("Kmer DB created");
            loadDB();
            kcount.loadHighCopyKmers();
            lg.verbose("Loaded high copy kmers");
            kcount.report(); // output
            kcount.cleanup(); // delete tmp files
            break;
            
        }
            
        case 2: { // union of multiple kmerdbs
            
            std::ifstream file;
            
            lg.verbose("Merging input databases");
            unsigned int numFiles = userInput.kmerDB.size(); // number of input kmerdbs
            
            uint8_t k = 0;
            
            for (uint16_t i = 0; i < numFiles; i++) {  // reads the kmer length from the index files checking consistency between kmerdbs
                
                file.open(userInput.kmerDB[i] + "/.index");
                std::string line;
                
                getline(file, line);
                file.close();
                
                k = stoi(line);
                
                if (k != stoi(line)) {
                    fprintf(stderr, "Cannot merge databases with different kmer length\n");
                    exit(1);
                }
                
            }
            
            if (k == 0) {
                fprintf(stderr, "Invalid kmer length (0)\n");
                exit(1);
            }
            
            userInput.kmerLen = k;
            kLen = userInput.kmerLen;
//            if (kLen < 31)
//                kPrefixLen = kLen;
            
            KDB kcount(userInput); // a new empty kmerdb with the specified kmer length
            
            lg.verbose("Kmer object generated. Merging.");
            kcount.kunion(); // union set
            kcount.report(); // output
            break;
        }
            
        default:
            fprintf(stderr, "Invalid mode\n");
            exit(1);
    }
    
}
