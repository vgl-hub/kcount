#ifndef KTREE_H
#define KTREE_H

class Knode {
    
    unsigned char* letter = NULL;
    Knode* children[4] = {NULL};

public:
    
    Knode(){};
    Knode(unsigned char* letter) : letter(letter) {};
    
    void set(unsigned char* c);
    
    friend class Ktree;
    
};

class Ktree {
    
    Knode* knodeRoot = NULL;
    unsigned short int KtreeH = 0;
    
    unsigned int totKmers = 0;
    unsigned int totKmersUnique = 0;
    
    Knode* nodes = new Knode[1000000000];
    unsigned long long int nodeCounter = 0;
    
    const unsigned char ctoi[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };
    
public:
    
    Ktree(InSequences& inSequences, unsigned short int k);
    
    ~Ktree();
    
    void delKnodeRecurse(Knode* current);
    
    void print2D(Knode* current, int space);
    
    void printKtree(Knode* root);
    
    void addChild(Knode* current, unsigned long long int pos, unsigned short int height, char* c);
    
    void addKmer(unsigned char* c);
    
};

#endif /* KTREE_H */


