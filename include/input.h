#ifndef INPUT_H
#define INPUT_H

struct UserInputKcount : UserInput {

    unsigned short int kmerLen = 21;

};

class Input {
    
    UserInputKcount userInput;
    
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
    
public:
    
    void load(UserInputKcount userInput);
    
    void read(InSequences& inSequence);
    
};

#endif /* INPUT_H */
