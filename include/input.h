#ifndef INPUT_H
#define INPUT_H

struct UserInputKcount : UserInput {

    unsigned short int kmerLen = 21; // default
    std::string outFile;

};

class Input {
    
    UserInputKcount userInput;
    
public:
    
    void load(UserInputKcount userInput);
    
    void read(bool mode);
    
};

#endif /* INPUT_H */
