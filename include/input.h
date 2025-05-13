#ifndef INPUT_H
#define INPUT_H

struct UserInputKcount : UserInput {
	uint32_t hashThreads = 7, writeThreads = 4; // small prime number, max
	int keepTmp = 0;
	uint8_t decompression_threads = 4, compression_threads = 6;
};

class Input {
    
    UserInputKcount userInput;
    
public:
    
    void load(UserInputKcount userInput);
    
    void loadDB();
    
    void read();
    
};

#endif /* INPUT_H */
