#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"
#include "input.h"

#include "kcount.h"

void Kcount::load(UserInputKcount& userInput) {
    
    InSequences inSequences;
    
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
                    
                    hashSegments(inSequences.getInSegments());
                    
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

                            appendReads(readBatch);

                            readBatch = new Sequences;

                        }

                        lg.verbose("Individual fastq sequence read: " + seqHeader);

                    }
                    
                    readBatch->batchN = seqPos/batchSize + 1;
                        
                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                    appendReads(readBatch);
                    
                    jobWait(threadPool);
                    
                    std::vector<Log> logs = inSequences.getLogs();
                    
                    //consolidate log
                    for (auto it = logs.begin(); it != logs.end(); it++) {
                        
                        it->print();
                        logs.erase(it--);
                        if(verbose_flag) {std::cerr<<"\n";};
                        
                    }

                    break;

                }
                    
            }
            
        }
        
    }
        
    count();
    
}

void Kcount::report(UserInputKcount& userInput) {
   
    std::cout<<"Total: "<<totKmers<<"\n";
    std::cout<<"Unique: "<<totKmersUnique<<"\n";
    std::cout<<"Distinct: "<<totKmersDistinct<<"\n";
    uint64_t missing = pow(4,k)-totKmersDistinct;
    std::cout<<"Missing: "<<missing<<"\n\n";
    
    const static phmap::flat_hash_map<std::string,int> string_to_case{
        {"hist",1}

    };
    
    // variable to handle output path and extension
    std::string path = rmFileExt(userInput.outSequence);
    std::string ext = getFileExt("." + userInput.outSequence);
    
    lg.verbose("Writing ouput: " + ext);
    
    // here we create a smart pointer to handle any kind of output stream
    std::unique_ptr<std::ostream> ostream;
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) { // this switch allows us to generate the output according to the input request and the unordered map. If the requested output format is not in the map we fall back to the undefined case (0)
            
        case 1: { // .hist
            
            std::ofstream ofs(userInput.outSequence);
            
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            
            printHist(ostream);
            
            break;
            
        }
        default: {
            
            ostream = std::make_unique<std::ostream>(std::cout.rdbuf());
            
            printHist(ostream);
            
        }
            
    }
    
}

bool Kcount::traverseInReads(Sequences* readBatch) { // traverse the read

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

void Kcount::appendReads(Sequences* readBatch) { // read a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch); });
    
    std::unique_lock<std::mutex> lck(mtx);
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
}

inline uint64_t Kcount::hash(uint8_t *kmer) {
    
    uint64_t fw = 0, rv = 0;
    
    for(uint8_t c = 0; c<k; ++c)
        fw += *kmer++ * pows[c];
    
    --kmer;
    
    for(uint8_t c = 0; c<k; ++c)
        rv += (3-(*kmer--)) * pows[c];
    
    return fw < rv ? fw : rv;
}

bool Kcount::countBuff(uint16_t m) {

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    buf64* thisBuf;
    
    phmap::flat_hash_map<uint64_t, uint64_t>* thisMap;
    
    for(buf64* buf : buffers) {
        
        thisBuf = &buf[m];
        
        thisMap = &map[m];
        
        uint64_t len = thisBuf->pos;
        
        for (uint64_t c = 0; c<len; ++c)
            ++(*thisMap)[thisBuf->seq[c]];
        
        delete[] thisBuf->seq;
        
    }
    
    return true;

}

bool Kcount::histogram(phmap::flat_hash_map<uint64_t, uint64_t>& map) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0;
    
    phmap::flat_hash_map<uint64_t, uint64_t> hist;
    
    for (auto pair : map) {
        
        if (pair.second == 1)
            ++kmersUnique;
        
        ++kmersDistinct;
        
        ++hist[pair.second];
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    totKmersUnique += kmersUnique;
    
    totKmersDistinct += kmersDistinct;
    
    for (auto pair : hist) {
        
        histogram1[pair.first] += pair.second;
        
    }
    
    return true;

}

void Kcount::printHist(std::unique_ptr<std::ostream>& ostream) {
    
    std::vector<std::pair<uint64_t, uint64_t>> table(histogram1.begin(), histogram1.end());
    std::sort(table.begin(), table.end());
    
    for (auto pair : table)
        *ostream<<pair.first<<"\t"<<pair.second<<"\n";

}

void Kcount::hashSequences(Sequences* readBatch) {
    
    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    buf64* buf = new buf64[mapCount];
    
    uint64_t kmers = 0;
    
    for (Sequence* sequence : readBatch->sequences) {
        
        uint64_t len = sequence->sequence->size(), kcount = len-k+1;
        
        if (len<k)
            continue;
        
        kmers += kcount;
        
        unsigned char* first = (unsigned char*)sequence->sequence->c_str();
        
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        uint64_t value, i, newSize;
        buf64* b;
        uint64_t* bufNew;
        
        for (uint64_t c = 0; c<kcount; ++c){
            
            value = hash(str+c);
            
            i = value / moduloMap;
            
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new uint64_t[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            b->seq[b->pos++] = value;
                        
        }
        
        delete[] str;
        
        threadLog.add("Processed sequence: " + sequence->header);
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    totKmers += kmers;
    
    buffers.push_back(buf);
    
    logs.push_back(threadLog);
    
}

void Kcount::hashSegments(std::vector<InSegment*>* segments) {
    
    buf64* buf = new buf64[mapCount];
    
    for (InSegment* segment : *segments) {
        
        uint64_t len = segment->getSegmentLen(), kcount = len-k+1;
        
        if (len<k)
            continue;
        
        totKmers += kcount;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        uint64_t value, i, newSize;
        buf64* b;
        uint64_t* bufNew;
        
        for (uint64_t c = 0; c<kcount; ++c){
            
            value = hash(str+c);
            
            i = value / moduloMap;
            
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new uint64_t[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            b->seq[b->pos++] = value;
                        
        }
        
        delete[] str;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
}

void Kcount::count() {
    
    lg.verbose("Counting with " + std::to_string(mapCount) + " maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return countBuff(m); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    lg.verbose("Generate histogram");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return histogram(map[m]); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
}
