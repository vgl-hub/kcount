#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cstdint>
#include <vector>
#include <thread>

#include "global.h"
#include "log.h"
#include "uid-generator.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"
#include "fastx.h"

std::string version = "0.0.1";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int hc_flag;
int hc_cutoff;
short int tabular_flag;
int cmd_flag;
int verbose_flag;
int outBubbles_flag;
int stats_flag;
int discoverPaths_flag;
int sortAlignment_flag;
int terminalAlignments_flag;
int maxThreads = 0;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;
std::vector<Log> logs;

struct KmerCounter {
	
	uint32_t k = 21, unique = 0, total = 0;
	std::unordered_map<std::string,uint32_t> kmers;
	
	KmerCounter(UserInput userInput) : k(userInput.kLen) {}
	
	void traverseInReads(Sequences *readBatch) {
		
		std::vector<Sequence*> sequences = readBatch->sequences;
		
		for (Sequence* sequence : sequences) {
			//std::cout<<"processing sequence: "<<sequence->header<<std::endl;
			if (sequence->sequence->size() >= k)
				appendSequence(sequence);
		}
		
		for (auto pair : kmers) {
			if (pair.second == 1)
				++unique;
		}
	}
	
	void consolidate() {}
	
	void appendSequence(Sequence* sequence) {
		for (uint32_t i = 0; i < sequence->sequence->size()-k+1; ++i) {
			std::string fw = sequence->sequence->substr(i, k);
			std::string rc = revCom(sequence->sequence->substr(i, k));
			std::string canonical = fw < rc ? fw : rc;
			std::cout<<canonical<<std::endl;
			++kmers[canonical];
			++total;
		}
	}
	
	void print() {
		std::cout<<"total kmers: "<<total<<std::endl;
		std::cout<<"unique kmers: "<<unique<<std::endl;
		std::cout<<"distinct kmers: "<<kmers.size()<<std::endl;
	}
};

int main(int argc, char **argv)
{
	
	UserInput userInput;
	short int c; // optarg
	bool isPipe = false; // to check if input is from pipe
	
	static struct option long_options[] = { // struct mapping long options
		{"input-sequences", required_argument, 0, 'r'},
		{"kmer-length", required_argument, 0, 'k'},		
		{0, 0, 0, 0}
	};
	
	while (true) { // loop through argv
		
		int option_index = 1;
		
		c = getopt_long(argc, argv, "-:k:r:",
						long_options, &option_index);
		
		if (c == -1) // exit the loop if run out of options
			break;
		switch (c) {
			case 'k': // kmer length
				
				if (!isNumber(optarg)) {
					fprintf(stderr, "input '%s' to option -%c must be a number\n", optarg, optopt);
					return EXIT_FAILURE;
				}
				
				userInput.kLen = atol(optarg);
				break;
			case 'r': // input reads
				if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
					userInput.pipeType = 'r'; // pipe input is a sequence
				}else{ // input is a regular file
					optind--;
					for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
						
						ifFileExists(argv[optind]);
						userInput.inFiles.push_back(argv[optind]);
					}
				}
				break;
		}
	}
	KmerCounter kmerCounter(userInput); // a new empty kmerdb with the specified kmer length
	unsigned int numFiles = userInput.inFiles.size();
	for (unsigned int i = 0; i < numFiles; ++i) // for each input read file
		loadSequences(userInput, kmerCounter, 'r', i); // specialized function to process reads into kmers as hashes
	
	kmerCounter.print();

	return 0;
}
