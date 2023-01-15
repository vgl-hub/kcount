#ifndef KCOUNT_H
#define KCOUNT_H

struct buf64 {
    uint64_t pos = 0, size = 100000;
    uint64_t *seq = new uint64_t[size];
};

class Kcount {
    
    std::vector<Log> logs;
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
    
    unsigned int batchSize = 10000;

    uint8_t k;
    
    uint64_t totKmers = 0, totKmersUnique = 0, totKmersDistinct = 0;
    
    const uint64_t mapCount = k < 28 ? pow(4,k/4) : pow(4,6);
    
    const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount;
    
    uint64_t* pows = new uint64_t[k];

    std::vector<buf64*> buffers;
    
    phmap::flat_hash_map<uint64_t, uint64_t>* map = new phmap::flat_hash_map<uint64_t, uint64_t>[mapCount];
    
    phmap::flat_hash_map<uint64_t, uint64_t> histogram1, histogram2;
    
    const uint8_t ctoi[256] = {
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    
public:
    
    Kcount(uint8_t k) : k(k) {
        
        for(uint8_t p = 0; p<k; ++p)
            pows[p] = (uint64_t) pow(4,p);
        
    };
    
    ~Kcount(){
        
        delete[] map;
        delete[] pows;
        
    }
    
    void load(UserInputKcount& userInput);
    
    bool traverseInReads(Sequences* readBatch);
    
    void appendReads(Sequences* readBatch);
    
    inline uint64_t hash(uint8_t* string);
    
    void hashSequences(Sequences* readBatch);
    
    void hashSegments(std::vector<InSegment*>* segments);
    
    void count();
    
    bool countBuff(uint16_t m);
    
    bool histogram(phmap::flat_hash_map<uint64_t, uint64_t>& map);
    
    void resizeBuff(buf64* buff);
    
    void printHist(std::unique_ptr<std::ostream>& ostream);
    
    void report(UserInputKcount& userInput);

};

#endif /* KCOUNT_H */
