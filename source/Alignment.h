#include <string>
#include <vector>
#ifndef KMER_H
#include "Kmer.h"
#endif // !KMER_H
//namespace that provides alignment functionality
namespace Alignment{
	//aligns reference and sequence k-mer, returns their local alignment and offset from reference
	std::pair<int,std::vector<std::string>> Align(const std::string &reference_kmer, const std::string &sequence_kmer);

	//backtracks in matrix vector based on operations
	void Backtrack(std::vector<int>operation, std::vector<int>matrix, int index,std::string const &ref, std::string const &seq, std::string &ref_align, std::string &seq_align);

	//initalizes matrix and operation vector
	void init_vectors(int k);

}