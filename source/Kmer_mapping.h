#include <string>
#include <vector>

class Kmer_mapping { 

public:
	int sequence_position;
	std::vector<int> reference_positions;

	Kmer_mapping(int x,int y);
	void add_position(int y);

};