#ifndef INPUT_H
#define INPUT_H

class Input {
    
    UserInput userInput;
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
    
public:
    
    void load(UserInput userInput);
    
    void read(InSequences& inSequence);
    
};

#endif /* INPUT_H */
