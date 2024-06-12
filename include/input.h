#ifndef INPUT_H
#define INPUT_H

struct UserInputKmap : UserInput {

};

class Input {
    
    UserInputKmap userInput;
    
public:
    
    void load(UserInputKmap userInput);
    
    void loadDB();
    
    void read(short unsigned int mode);
    
};

#endif /* INPUT_H */
