#ifndef INPUT_H
#define INPUT_H

struct UserInputKmap : UserInput {

    unsigned short int kmerLen = 21; // default
    std::string outFile;

};

class Input {
    
    UserInputKmap userInput;
    
public:
    
    void load(UserInputKmap userInput);
    
    void read(bool mode);
    
};

#endif /* INPUT_H */
