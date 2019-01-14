#include "Alignment.h"
#include <set>
#include <unordered_map>

class Kmer;

namespace mapping {
	//maps sequence to reference, returns position if successful, else -99999
	int map_to_reference(std::vector<Kmer>* reference, std::vector<Kmer>* sequence);
	//maps sequence to reference using map of odrering number and position, returns position if successful, else -99999
	int alternative_mapping(const std::unordered_map<int,int> &reference,const std::vector<Kmer> & sequence);
	//checks if char count difference is less than 4
	bool check_match(const std::string &k1, const std::string &k2);
};


