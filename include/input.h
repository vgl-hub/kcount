#ifndef INPUT_H
#define INPUT_H

struct UserInputKcount : UserInput {
	uint32_t hashThreads = 7, writeThreads = 1; // small prime number
	int keepTmp = 0;
};

class Input {
    
    UserInputKcount userInput;
    
public:
    
    void load(UserInputKcount userInput);
    
    void loadDB();
    
    void read();
    
};

#endif /* INPUT_H */
