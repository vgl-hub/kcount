#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"

#include "ktree.h"

Knode* Knode::contains(char c){
    
    switch(c) {
        case 'A':
            if (this->A != NULL)
                return this->A;
            break;
        case 'C':
            if (this->C != NULL)
                return this->C;
            break;
        case 'G':
            if (this->G != NULL)
                return this->G;
            break;
        case 'T':
            if (this->T != NULL)
                return this->T;
            break;
    }
    
    return NULL;
    
}

void Ktree::print2D(Knode* current, int space) {
    
    if (current == NULL)
        return;
    
    space += 3;
    
    // Process right children first
    print2D(current->A, space);
    print2D(current->C, space);
    
    printf("\n%*s%c\n", space-3, "", *(current->letter));
    
    // Process left children
    print2D(current->G, space);
    print2D(current->T, space);
    
    return;
    
}

void Ktree::printKtree(Knode* root) {
    print2D(root, 0);
}

void Ktree::addChild(Knode* current, unsigned long long int pos, unsigned short int height, char* c) {
    
    if (height>=KtreeH)
        return;
    
    switch (*(c+pos)) {
            
        case 'A':
            
            if (current->A == NULL)
                current->A = new Knode(height, c+pos);
            addChild(current->A, ++pos, ++height, c);
            break;
            
        case 'C':
            
            if (current->C == NULL)
                current->C = new Knode(height, c+pos);
            addChild(current->C, ++pos, ++height, c);
            break;
            
        case 'G':
            
            if (current->G == NULL)
                current->G = new Knode(height, c+pos);
            addChild(current->G, ++pos, ++height, c);
            break;
            
        case 'T':
            
            if (current->T == NULL)
                current->T = new Knode(height, c+pos);
            addChild(current->T, ++pos, ++height, c);
            
    }
    
    return;
    
}

void Ktree::addKmer(char* c) {
    addChild(knodeRoot, 0, 0, c);
}

Ktree::Ktree(InSequences& inSequences, unsigned short int k) {
    
    KtreeH = k;
    
    lg.verbose("Started ktree construction");
    
    knodeRoot = new Knode(0, new char('0'));
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    
    for (InSegment* segment : *segments) {
        
        long long int len = segment->getSegmentLen()-k+1;
        
        char* first = segment->first();
        
        for (long long int c = 0; c<len; ++c) {
            
            addKmer(first+c);
            
        }
        
    }
            
    printKtree(knodeRoot);
    
}

void Ktree::delKnodeRecurse(Knode* current) {
    
    if (current == NULL)
        return;

    delKnodeRecurse(current->T);
    delKnodeRecurse(current->G);
    delKnodeRecurse(current->C);
    delKnodeRecurse(current->A);
    
    delete current;
    
    return;
    
}

Ktree::~Ktree() {
    
    delete knodeRoot->letter;
    delKnodeRecurse(knodeRoot);
    
}
