#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"

#include "kcount.h"

bool Kcount::traverseInReads(Sequences* readBatch) { // traverse the read

    Log threadLog;
    
    threadLog.setId(readBatch->batchN);

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

void Kcount::appendReads(Sequences* readBatch) { // read a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch); });
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    lck.unlock();
    
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

void Kcount::printHist() {
    
    std::vector<std::pair<uint64_t, uint64_t>> table(histogram1.begin(), histogram1.end());
    std::sort(table.begin(), table.end());
    
    for (auto pair : table)
        std::cout<<pair.first<<"\t"<<pair.second<<"\n";

}

void Kcount::hashSequences(Sequences* readBatch) {
    
    buf64* buf = new buf64[mapCount];
    
    for (Sequence* sequence : readBatch->sequences) {
        
        uint64_t len = sequence->sequence->size(), kcount = len-k+1;
        
        if (len<k)
            continue;
        
        totKmers += kcount;
        
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
        
        lg.verbose("Processed sequence: " + sequence->header);
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
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
    
    std::cout<<"Total: "<<totKmers<<"\n";
    std::cout<<"Unique: "<<totKmersUnique<<"\n";
    std::cout<<"Distinct: "<<totKmersDistinct<<"\n";
    uint64_t missing = pow(4,k)-totKmersDistinct;
    std::cout<<"Missing: "<<missing<<"\n\n";
    
    printHist();
    
}
