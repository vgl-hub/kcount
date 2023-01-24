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

void Input::read(bool mode) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    if (mode == 0) {
        
        Kmap<uint64_t> kcount(userInput.kmerLen);
        
        lg.verbose("Kmer object generated");
        
        kcount.convert(userInput);
        
        kcount.count();
        
        kcount.report(userInput);
        
    }else{
        
        std::ifstream file;

        file.open(userInput.iSeqFileArg + "/.index");
        std::string line;
        
        getline(file, line);
        file.close();
        
        short unsigned int k = stoi(line);
        
        Kmap<uint64_t> kcount(k);
        
        lg.verbose("Kmer object generated");
        
        kcount.load(userInput);
        
        kcount.report(userInput);
        
    }
    
}
