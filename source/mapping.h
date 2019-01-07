#include "Alignment.h"
#include <set>
#include <unordered_map>

class Kmer;

namespace mapping {

	int map_to_reference(std::vector<Kmer>* reference, std::vector<Kmer>* sequence);
	int alternative_mapping(const std::unordered_map<int,int> &reference,const std::vector<Kmer> & sequence);
	bool check_match(const std::string &k1, const std::string &k2);
};


