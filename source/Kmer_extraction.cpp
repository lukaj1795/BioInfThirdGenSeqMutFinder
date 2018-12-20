#include"Kmer_extraction.h"
#include"Genome.h"
#include"Kmer.h"
#include <algorithm>
#include <vector>

bool sort_Kmers(Kmer k1, Kmer k2)
{
	return (k1.ordering_number_for_string < k2.ordering_number_for_string);
}
std::vector<Kmer> Kmer_extraction::extract(Genome *sequence) {
	
	int n = sequence->genomeString.length();
	std::vector<Kmer> all_kmers;

	/*the length of String*/
	for (int j = 0; j <= n - window_size; ++j) {
		std::vector<Kmer> minimizer_lookup;
		/* cycle over the window till k-mers of length, k, can still be made */
		for (int i = j; i < j+ w; i = i + 1) {
			minimizer_lookup.push_back(Kmer(sequence->genomeString.substr(i, k), i, sequence->identifier));
		}
		std::sort(minimizer_lookup.begin(), minimizer_lookup.end(), sort_Kmers);
		bool found = false;
		// Iterate over all elements in Vector
		for (auto & elem : all_kmers){
			if (elem.string == minimizer_lookup.front().string)
			{
				found = true;
				break;
			}
		}
		if (!found) {
			all_kmers.push_back(minimizer_lookup.front());
		}
	}
	return all_kmers;
};


Kmer_extraction::Kmer_extraction(int x, int y) {
	w = x;
	k = y;
	window_size = w + k - 1;
}
