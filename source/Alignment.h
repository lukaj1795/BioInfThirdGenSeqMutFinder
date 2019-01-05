#include <string>
#include <vector>
#ifndef KMER_H
#include "Kmer.h"
#endif // !KMER_H

namespace Alignment{

	std::pair<int,std::vector<std::string>> Align(Kmer reference, Kmer sequences);

	void Backtrack(std::vector<int>operation, std::vector<int>matrix, int index,std::string const &ref, std::string const &seq, std::string &ref_align, std::string &seq_align);

	void init_vectors(int k);

}