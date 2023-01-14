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

#include "kcount.h"
#include "input.h"

void Input::load(UserInputKcount userInput) {
    
    this->userInput = userInput;
    
}


void Input::read(InSequences& inSequences) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    Kcount kcount(userInput.kmerLen);
    
    threadPool.init(maxThreads); // initialize threadpool
    
    stream = streamObj.openStream(userInput, 'f');
    
    if (!userInput.iSeqFileArg.empty() || userInput.pipeType == 'f') {
        
        StreamObj streamObj;
        
        stream = streamObj.openStream(userInput, 'f');
        
        if (stream) {
            
            switch (stream->peek()) {
                    
                case '>': {
                    
                    stream->get();
                    
                    while (getline(*stream, newLine)) {
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
                        std::string* inSequence = new std::string;
                        
                        getline(*stream, *inSequence, '>');
                        
                        lg.verbose("Individual fasta sequence read");
                        
                        Sequence* sequence = new Sequence{seqHeader, seqComment, inSequence, NULL};
                        
                        sequence->seqPos = seqPos; // remember the order
                        
                        inSequences.appendSequence(sequence);
                        
                        seqPos++;
                        
                    }
                    
                    jobWait(threadPool);
                    
                    if(verbose_flag) {std::cerr<<"\n\n";};
                    
                    std::vector<Log> logs = inSequences.getLogs();
                    
                    //consolidate log
                    for (auto it = logs.begin(); it != logs.end(); it++) {
                        
                        it->print();
                        logs.erase(it--);
                        if(verbose_flag) {std::cerr<<"\n";};
                        
                    }
                    
                    kcount.hashSegments(inSequences.getInSegments());
                    
                    break;
                    
                }
                    
                case '@': {
                    
                    Sequences* readBatch = new Sequences;

                    while (getline(*stream, newLine)) { // file input

                        newLine.erase(0, 1);

                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }

                        std::string* inSequence = new std::string;
                        getline(*stream, *inSequence);

                        getline(*stream, newLine);
                        
                        ignore(*stream, '\n');

                        readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence});
                        seqPos++;

                        if (seqPos % batchSize == 0) {

                            readBatch->batchN = seqPos/batchSize;
                            
                            lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                            kcount.appendReads(readBatch);

                            readBatch = new Sequences;

                        }

                        lg.verbose("Individual fastq sequence read: " + seqHeader);

                    }
                    
                    readBatch->batchN = seqPos/batchSize + 1;
                        
                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                    kcount.appendReads(readBatch);

                    break;

                }
                    
            }
            
        }
        
    }
    
    jobWait(threadPool);
        
    kcount.count();
	
	threadPool.join();
    
}
