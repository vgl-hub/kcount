#ifndef INPUT_H
#define INPUT_H

struct UserInputKcount : UserInput {

};

class Input {
    
    UserInputKcount userInput;
    
public:
    
    void load(UserInputKcount userInput);
    
    void loadDB();
    
    void read(short unsigned int mode);
    
};

#endif /* INPUT_H */
