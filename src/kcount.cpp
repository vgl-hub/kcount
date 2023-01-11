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

inline uint64_t Kcount::hash(uint8_t *kmer) {
    
    uint64_t fw = 0, rv = 0;
    
    for(uint8_t c = 0; c<k; ++c)
        fw += *kmer++ * pows[c];
    
    --kmer;
    
    for(uint8_t c = 0; c<k; ++c)
        rv += (3-(*kmer--)) * pows[c];
    
    return fw < rv ? fw : rv;
}

bool Kcount::countBuff(buf64* buf, phmap::flat_hash_map<uint64_t, uint64_t>& map) {

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    uint64_t len = buf->pos;
    
    for (uint64_t c = 0; c<len; ++c)
        ++map[buf->seq[c]];
    
    return true;

}

bool Kcount::countUnique(phmap::flat_hash_map<uint64_t, uint64_t>& map) {
    
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

void Kcount::count(std::vector<InSegment*>* segments) {
    
    uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount;
    
    lg.verbose("Counting with " + std::to_string(mapCount) + " maps");
    
    for (InSegment* segment : *segments) {
        
        if (segment->getSegmentLen()<k)
            continue;
        
        uint64_t len = segment->getSegmentLen(), kcount = len-k+1;
        
        totKmers += kcount;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        uint8_t* str = new uint8_t[segment->getSegmentLen()];
        
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
    
    lg.verbose("Populating maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return countBuff(&buf[m], map[m]); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    lg.verbose("Counting unique kmers");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return countUnique(map[m]); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    std::cout<<"Total: "<<totKmers<<"\n";
    std::cout<<"Unique: "<<totKmersUnique<<"\n";
    std::cout<<"Distinct: "<<totKmersDistinct<<"\n";
    uint64_t missing = pow(4,k)-totKmersDistinct;
    std::cout<<"Missing: "<<missing<<"\n\n";
    
    printHist();
    
}
