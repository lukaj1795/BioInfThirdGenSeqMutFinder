#include "Alignment.h"
#include <set>

class Kmer;

namespace mapping {

	int map_to_reference(std::vector<Kmer>* reference, std::vector<Kmer>* sequence);
	//int alternative_mapping(std::vector<Kmer>* reference, std::vector<Kmer>* sequence);
	bool check_match(const Kmer &k1, const Kmer &k2);
};


