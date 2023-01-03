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

void Ktree::print2D(Knode* current, int space) {
    
    if (current == NULL)
        return;
    
    space += 3;
    
    // Process right children first
    print2D(current->children[0], space);
    print2D(current->children[1], space);
    
    printf("\n%*s%c\n", space-3, "", *(current->letter));
    
    // Process left children
    print2D(current->children[2], space);
    print2D(current->children[3], space);
    
    return;
    
}

void Ktree::printKtree(Knode* root) {

    print2D(root, 0);
    
}

//void Ktree::addChild(Knode* current, unsigned long long int pos, unsigned short int height, char* c) {
//
//    if (height>=KtreeH)
//        return;
//
//    switch (*(c+pos)) {
//
//        case 'A':
//
//            if (current->children[0] == NULL)
//                current-> = new Knode(c+pos);
//            addChild(current->A, ++pos, ++height, c);
//            break;
//
//        case 'C':
//
//            if (current->C == NULL)
//                current->C = new Knode(c+pos);
//            addChild(current->C, ++pos, ++height, c);
//            break;
//
//        case 'G':
//
//            if (current->G == NULL)
//                current->G = new Knode(c+pos);
//            addChild(current->G, ++pos, ++height, c);
//            break;
//
//        case 'T':
//
//            if (current->T == NULL)
//                current->T = new Knode(c+pos);
//            addChild(current->T, ++pos, ++height, c);
//
//    }
//
//    return;
//
//}

//void Ktree::addKmer(char* c) {
//    addChild(knodeRoot, 0, 0, c);
//}

//void Ktree::addKmer(char* c) {
//
//    Knode* current = knodeRoot;
//
//    for (unsigned short int height = 0; height<KtreeH; height++) {
//
//        current = current->alphabet[*(c+height)];
//
//        if (current == NULL)
//            current = new Knode(c+height);
//        printf("%c\n", *(current->letter));
//
//    }
//
//}

void Knode::set(unsigned char* c){
    
    this->letter = c;
    
}

void Ktree::addKmer(unsigned char* c) {
    
    Knode* current = knodeRoot;
    
    for (unsigned short int height = 0; height<KtreeH; height++) {

        if (current->children[ctoi[*(c+height)]] == NULL) {
            
            current->children[ctoi[*(c+height)]] = &nodes[nodeCounter++];
            
            current->set(c+height);
            
            if (height+1==KtreeH)
                ++totKmersUnique;
            
        }
        current = current->children[ctoi[*(c+height)]];
        
    }
    
    ++totKmers;
    
}

Ktree::Ktree(InSequences& inSequences, unsigned short int k) {
    
    KtreeH = k;
    
    lg.verbose("Started ktree construction");
    
    knodeRoot = new Knode(new unsigned char('0'));
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    
    for (InSegment* segment : *segments) {
        
        lg.verbose("Processing segment: " + segment->getSeqHeader());
        
        if (segment->getSegmentLen()<k) {
            lg.verbose("Segment shorted thank k. skipping");
            continue;
        }
        
        unsigned long long int len = segment->getSegmentLen()-k+1;

        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        for (unsigned long long int c = 0; c<len; ++c) {

            addKmer(first+c);
            
        }
        
    }
            
//    printKtree(knodeRoot);
    
    std::cout<<"Total kmers: "<<totKmers<<std::endl;
    std::cout<<"Unique kmers: "<<totKmersUnique<<std::endl;
    
}

void Ktree::delKnodeRecurse(Knode* current) {
    
    if (current == NULL)
        return;

    delKnodeRecurse(current->children[3]);
    delKnodeRecurse(current->children[2]);
    delKnodeRecurse(current->children[1]);
    delKnodeRecurse(current->children[0]);
    
    delete current;
    
    return;
    
}

Ktree::~Ktree() {
    
    delete knodeRoot->letter;
    delete[] nodes;
//    delKnodeRecurse(knodeRoot);
    
}
