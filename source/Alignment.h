#include <string>
#include <vector>

#ifndef KMER_H
#include "Kmer.h"
#endif // !KMER_H

namespace Alignment{

	std::string Align(std::vector<Kmer> reference, std::vector<Kmer> sequences);

	void Backtrack(int *operation, int *matrix, int column, int index,std::string const &ref, std::string const &seq, std::string &ref_align, std::string &seq_align);


}