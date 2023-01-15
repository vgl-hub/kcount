#ifndef INPUT_H
#define INPUT_H

struct UserInputKcount : UserInput {

    unsigned short int kmerLen = 21; // default

};

class Input {
    
    UserInputKcount userInput;
    
public:
    
    void load(UserInputKcount userInput);
    
    void read();
    
};

#endif /* INPUT_H */
