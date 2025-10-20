#ifndef KCOUNT_H
#define KCOUNT_H

#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/bgzf.h>

class KDB : public Kmap<KDB, UserInputKcount, Key, uint8_t, uint32_t> { // CRTP
    UserInputKcount userInput;
public:
    KDB(UserInputKcount& userInput) : Kmap{userInput}, userInput(userInput) {
        DBextension = "kc";
    }
};

#define BUFFER_RESERVE_SIZE 2097152 // 2MB preallocation for efficiency
template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::initBuffering() {
	
	Init_Genes_Package(k, sLen);
	std::string kmerBatch;
	kmerBatch.reserve(BUFFER_RESERVE_SIZE); // Preallocate memory for efficiency
	Distribution_Bundle* bundle = Begin_Distribution(bufferFiles);
	
	uint32_t batchSize = 10000000; // number of bases processed by a thread
	htsThreadPool tpool_read; // htslib threadpool pointer
	std::string newLine, seqHeader, seqComment, line, bedHeader;
	std::size_t numFiles = userInput.inFiles.size();
	uint32_t batchN = 0;
	uint64_t processedLength = 0;
	lg.verbose("Processing " + std::to_string(numFiles) + " files");
	
	const static phmap::flat_hash_map<std::string,int> string_to_case{
		{"fasta",1},
		{"fa",1},
		{"fasta.gz",1},
		{"fa.gz",1},
		{"fastq",1},
		{"fq",1},
		{"fastq.gz",1},
		{"fq.gz",1},
		{"bam",1},
		{"cram",1}
	};
	
	for (uint32_t i = 0; i < numFiles; i++) {
		
		std::string file = userInput.file('r', i);
		std::string ext = getFileExt(file);
		lg.verbose("Processing file: " + file);
		
		switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
				
			case 1: { // fa*[.gz], bam, cram
				
				samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); // open file
				bam_hdr_t *bamHdr = sam_hdr_read(fp_in); // read header
				bam1_t *bamdata = bam_init1(); // initialize an alignment
				std::string inSequence; // new sequence that we can reuse
				int64_t pos_before = bgzf_tell(fp_in->fp.bgzf) >> 16;

				int64_t compressed_bytes_read = 0, totalKmers = 0;
				
				tpool_read = {NULL, 0};
				tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
				if (tpool_read.pool)
					hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
				else
					lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");
				
				while(sam_read1(fp_in,bamHdr,bamdata) > 0) {
					
					int64_t pos_after = bgzf_tell(fp_in->fp.bgzf) >> 16; // keep track of the original file size to estimate compression
					compressed_bytes_read += pos_after - pos_before;
					pos_before = pos_after;
					
					uint32_t len = bamdata->core.l_qseq; // length of the read
					totalKmers += len;
					uint8_t *seq = bam_get_seq(bamdata); // seq string
					
					inSequence.clear();
					inSequence.resize(len);
					for(uint32_t i=0; i<len; ++i)
						inSequence.at(i) = seq_nt16_str[bam_seqi(seq,i)]; // gets nucleotide id and converts them into IUPAC id
					
					kmerBatch += inSequence;
					kmerBatch.push_back('N');
					processedLength += inSequence.size();
					
					if (processedLength > batchSize) {
						lg.verbose("Processing batch N: " + std::to_string(batchN++));
						lg.verbose("Found " + std::to_string(totalKmers) + " total kmers (extracted from " + std::to_string(compressed_bytes_read) + " bytes)");
						Distribute_Sequence(const_cast<char*>(kmerBatch.data()), kmerBatch.size(), bundle);
						processedLength = 0;
						kmerBatch.clear();
					}
				}
				lg.verbose("Processing batch N: " + std::to_string(batchN++));
				lg.verbose("Found " + std::to_string(totalKmers) + " total kmers (extracted from " + std::to_string(compressed_bytes_read) + " bytes)");
				Distribute_Sequence(const_cast<char*>(kmerBatch.data()), kmerBatch.size(), bundle);
				
				bam_destroy1(bamdata);
				sam_close(fp_in);
				if (tpool_read.pool)
					hts_tpool_destroy(tpool_read.pool);
				break;
			}
			default: {
				fprintf(stderr, "cannot recognize input (must be: fasta, fastq, bam, cram).\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	End_Distribution(bundle);
	finalize();
}


#endif /* KCOUNT_H */
