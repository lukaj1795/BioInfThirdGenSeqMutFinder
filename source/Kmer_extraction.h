#include <vector>

class Genome;
class Kmer;

class Kmer_extraction {

private:

	int w; /*w consecutive k-mers*/
	int k; /*length of k-mer, number of letters in one k-mer*/
	int window_size; /*a set of w consecutive k-mers covers a string=window of exactly w+k-1 letters*/

public:
	//extracts k-mers from sequence
	std::vector<Kmer> extract(Genome *sequence);

	Kmer_extraction(int x, int y);
	//extracts k-mers from sequence complement 
	std::vector<Kmer> extract_complement(Genome * sequence);

};
