#include <string>
#include <vector>

class Kmer_extraction {

public:

	int w;
	int k;
	std::vector<std::string> all_kmers;

	std::vector<std::string> extract(const std::string& sequence);


	Kmer_extraction(int x, int y) {
		w = x;
		k = y;
	}

};
