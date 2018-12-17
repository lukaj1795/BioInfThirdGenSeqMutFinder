#include"Kmer_extraction.h"


std::vector<std::string> Kmer_extraction::extract(const std::string& sequence) {

	int n = sequence.length();
	std::vector<std::string> all_kmers;
	/* cycle over the length of String till k-mers of length, k, can still be made */
	for (int i = 0; i < n - k + 1; i = i + w) {
		all_kmers.push_back(sequence.substr(i,k));
	}


	return all_kmers;

};

