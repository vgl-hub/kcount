#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <getopt.h>
#include <fstream>

#include "global.h"
#include "log.h"
#include "uid-generator.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"

#include "input.h"
#include "main.h"

std::string version = "0.1.0";

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

uint32_t kLen = 63;
uint8_t kPrefixLen = 31;
Buf<uint8_t> *seqBuf, *seqBuf2;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;
std::vector<Log> logs;

void printHelp() {
    
    printf("kcount [mode] -h\nfor additional help.\n");
    printf("\nModes:\n");
    printf("count\n");
    printf("load\n");
    printf("union\n");
    exit(0);
    
}

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    bool isPipe = false; // to check if input is from pipe
    std::string cmd;
    UserInputKcount userInput; // initialize input object
    
    if (argc == 1 || argc == 2) // kcount with no arguments
        printHelp();
    const static phmap::parallel_flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"count",0},
        {"load",1},
        {"union",2}
    };
    
    auto got = string_to_case.find(argv[1]);
    if (got != string_to_case.end()) {
        userInput.mode = got->second;
    }else{
        fprintf(stderr, "mode %s does not exist. Terminating\n", argv[1]);
        return EXIT_FAILURE;
    }
    
    switch (userInput.mode) {
        case 0: { // count
            
            static struct option long_options[] = { // struct mapping long options
                {"input-sequences", required_argument, 0, 'r'},
                {"kmer-length", required_argument, 0, 'k'},
                {"out-format", required_argument, 0, 'o'},
                
                {"threads", required_argument, 0, 'j'},
                {"verbose", no_argument, &verbose_flag, 1},
                {"cmd", no_argument, &cmd_flag, 1},
                
                {"help", no_argument, 0, 'h'},
                
                {0, 0, 0, 0}
            };
            
            while (true) { // loop through argv
                
                int option_index = 1;
                
                c = getopt_long(argc, argv, "-:k:j:o:r:h",
                                long_options, &option_index);
                
                if (c == -1) { // exit the loop if run out of options
                    break;
                    
                }
                
                switch (c) {
                    case ':': // handle options without arguments
                        switch (optopt) { // the command line option last matched
                            case 'b':
                                break;
                                
                            default:
                                fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                                return EXIT_FAILURE;
                        }
                        break;
                    default: // handle positional arguments
                        
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        
                        break;
                    case '?': // unrecognized option
                        if (optopt != 0)
                            fprintf(stderr, "Unrecognized option: %c\n", optopt);
                        else
                            fprintf(stderr, "Unrecognized option: %s\n", argv[optind-1]);
                        printHelp();
                        break;
                    case 'k': // kmer length
                        
                        if (!isNumber(optarg)) {
                            fprintf(stderr, "input '%s' to option -%c must be a number\n", optarg, optopt);
                            return EXIT_FAILURE;
                        }
                        
                        userInput.kmerLen = atoi(optarg);
                        kLen = userInput.kmerLen;
                        if (kLen < 31)
                            kPrefixLen = kLen;
                        break;
                        
                    case 'j': // max threads
                        maxThreads = atoi(optarg);
                        break;
                        
                    case 'o': // handle output (file or stdout)
                        userInput.outFile = optarg;
                        break;
                        
                    case 'r': // input reads
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'r'; // pipe input is a sequence
                        }else{ // input is a regular file
                            optind--;
                            for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
                                
                                ifFileExists(argv[optind]);
                                userInput.inReads.push_back(argv[optind]);
                            }
                        }
                        break;
                    case 'v': // software version
                        printf("kcount v%s\n", version.c_str());
                        printf("Giulio Formenti giulio.formenti@gmail.com\n");
                        exit(0);
                        
                    case 'h': // help
                        printf("kcount count [options]\n");
                        printf("\nOptions:\n");
                        printf("\t-r --input-sequences sequence input files (fasta,fastq).\n");
                        printf("\t-k --kmer-length length of kmers (default: 63).\n");
                        printf("\t-j --threads <n> numbers of threads (default: max).\n");
                        printf("\t-o --out-format generates various kinds of outputs (currently supported: .hist .kc).\n");
                        printf("\t--verbose verbose output.\n");
                        printf("\t--cmd print $0 to stdout.\n");
                        exit(0);
                }
                
                if    (argc == 2 || // handle various cases in which the output should include default outputs
                       (argc == 3 && pos_op == 2) ||
                       (argc == 4 && pos_op == 3)) {
                    
                }
                
            }
            
            if (userInput.inReads.size() == 0) {
                fprintf(stderr, "At least one input file required (-f).\n");
                return EXIT_FAILURE;
            }
            
            break;
        }
        case 1: { // load
            
            static struct option long_options[] = { // struct mapping long options
                {"database", required_argument, 0, 'd'},
                {"out-format", required_argument, 0, 'o'},
                
                {"threads", required_argument, 0, 'j'},
                {"verbose", no_argument, &verbose_flag, 1},
                {"cmd", no_argument, &cmd_flag, 1},
                {"help", no_argument, 0, 'h'},
                
                {0, 0, 0, 0}
            };
            
            while (true) { // loop through argv
                
                int option_index = 1;
                
                c = getopt_long(argc, argv, "-:d:j:o:h",
                                long_options, &option_index);
                
                if (c == -1) { // exit the loop if run out of options
                    break;
                    
                }
                
                switch (c) {
                    case ':': // handle options without arguments
                        switch (optopt) { // the command line option last matched
                            case 'b':
                                break;
                                
                            default:
                                fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                                return EXIT_FAILURE;
                        }
                        break;
                    default: // handle positional arguments
                        
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        
                        break;
                    case '?': // unrecognized option
                        if (optopt != 0)
                            fprintf(stderr, "Unrecognized option: %c\n", optopt);
                        else
                            fprintf(stderr, "Unrecognized option: %s\n", argv[optind-1]);
                        printHelp();
                        break;
                    case 'd': // input DBs
                        
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'f'; // pipe input is a sequence
                        }else{ // input is a regular file
                            
                            optind--;
                            for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
                                
                                ifFileExists(argv[optind]);
                                userInput.kmerDB.push_back(argv[optind]);
                            }
                        }
                        break;
                    case 'j': // max threads
                        maxThreads = atoi(optarg);
                        break;
                        
                    case 'o': // handle output (file or stdout)
                        userInput.outFile = optarg;
                        break;
                        
                    case 'h': // help
                        printf("kcount load [options]\n");
                        printf("\nOptions:\n");
                        printf("\t-d --database kmer database to load.\n");
                        printf("\t-j --threads <n> numbers of threads (default: max).\n");
                        printf("\t-o --out-format generates various kinds of outputs (currently supported: .hist .kc).\n");
                        printf("\t--verbose verbose output.\n");
                        printf("\t--cmd print $0 to stdout.\n");
                        exit(0);
                }
                
                if    (argc == 2 || // handle various cases in which the output should include default outputs
                       (argc == 3 && pos_op == 2) ||
                       (argc == 4 && pos_op == 3)) {
                    
                }
                
            }
            
            if (userInput.kmerDB.size() != 1) {
                fprintf(stderr, "One database required (-d).\n");
                return EXIT_FAILURE;
            }
            
            break;
        }
        case 2: { // union
            
            static struct option long_options[] = { // struct mapping long options
                {"databases", required_argument, 0, 'd'},
                {"out-format", required_argument, 0, 'o'},
                
                {"threads", required_argument, 0, 'j'},
                {"verbose", no_argument, &verbose_flag, 1},
                {"cmd", no_argument, &cmd_flag, 1},
                {"help", no_argument, 0, 'h'},
                
                {0, 0, 0, 0}
            };
            
            while (true) { // loop through argv
                
                int option_index = 1;
                
                c = getopt_long(argc, argv, "-:d:j:o:h",
                                long_options, &option_index);
                
                if (c == -1) { // exit the loop if run out of options
                    break;
                    
                }
                
                switch (c) {
                    case ':': // handle options without arguments
                        switch (optopt) { // the command line option last matched
                            case 'b':
                                break;
                                
                            default:
                                fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                                return EXIT_FAILURE;
                        }
                        break;
                    default: // handle positional arguments
                        
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        
                        break;
                    case '?': // unrecognized option
                        if (optopt != 0)
                            fprintf(stderr, "Unrecognized option: %c\n", optopt);
                        else
                            fprintf(stderr, "Unrecognized option: %s\n", argv[optind-1]);
                        printHelp();
                        break;
                    case 'd': // input DBs
                        
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'f'; // pipe input is a sequence
                        }else{ // input is a regular file
                            
                            optind--;
                            for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
                                
                                ifFileExists(argv[optind]);
                                userInput.kmerDB.push_back(argv[optind]);
                            }
                        }
                        break;
                    case 'j': // max threads
                        maxThreads = atoi(optarg);
                        break;
                        
                    case 'o': // handle output (file or stdout)
                        userInput.outFile = optarg;
                        break;
                        
                    case 'h': // help
                        printf("kcount union [options]\n");
                        printf("\nOptions:\n");
                        printf("\t-d --databases kmer database to merge.\n");
                        printf("\t-j --threads <n> numbers of threads (default: max).\n");
                        printf("\t-o --out-format generates various kinds of outputs (currently supported: .hist .kc).\n");
                        printf("\t--verbose verbose output.\n");
                        printf("\t--cmd print $0 to stdout.\n");
                        exit(0);
                }
                
                if    (argc == 2 || // handle various cases in which the output should include default outputs
                       (argc == 3 && pos_op == 2) ||
                       (argc == 4 && pos_op == 3)) {
                    
                }
                
            }
            
            if (userInput.kmerDB.size() < 2) {
                fprintf(stderr, "At least two databases required (-d).\n");
                return EXIT_FAILURE;
            }
            break;
        }
        default: {
            fprintf(stderr, "Unrecognized mode: %s\n", argv[1]);
            printHelp();
        }
    }
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
        
	Input in;
	in.load(userInput); // load user input
	lg.verbose("Loaded user input");
    maxMem = (userInput.maxMem == 0 ? get_mem_total(3) * 0.9 : userInput.maxMem); // set memory limit
    
    threadPool.init(maxThreads); // initialize threadpool
    
    in.read(); // read input
        
    threadPool.join(); // join threads

    exit(EXIT_SUCCESS);
	
}
