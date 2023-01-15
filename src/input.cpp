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
#include "kcount.h"

void Input::load(UserInputKcount userInput) {
    
    this->userInput = userInput;
    
}

void Input::read() {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    Kcount kcount(userInput.kmerLen);
    
    lg.verbose("Kmer object generated");
    
    kcount.load(userInput);

    kcount.report(userInput);
    
}
