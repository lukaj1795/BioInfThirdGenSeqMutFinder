#include"Kmer_extraction.h"
#include"Genome.h"
#include"Kmer.h"
#include <algorithm>
#include <vector>

bool sort_Kmers(Kmer k1, Kmer k2)
{
	return (k1.ordering_number_for_string < k2.ordering_number_for_string);
}
bool sort_pozicija(Kmer k1, Kmer k2)
{
	return (k1.position < k2.position);
}
std::vector<Kmer> Kmer_extraction::extract(Genome *sequence) {

	int n = sequence->genomeString.length();
	std::vector<Kmer> all_kmers;

	/*the length of String*/
	for (int j = 0; j <= n - window_size; ++j) {
		std::vector<Kmer> minimizer_lookup;
		/* cycle over the window till k-mers of length, k, can still be made */
		for (int i = j; i < j + w; i = i + 1) {
			minimizer_lookup.push_back(Kmer(sequence->genomeString.substr(i, k), i, sequence->identifier));
		}
		/*only smallest k-mer will become minimizer of the given window*/
		std::sort(minimizer_lookup.begin(), minimizer_lookup.end(), sort_Kmers);
		bool found = false;
		// Iterate over all elements in Vector and search if the k-mer is already contained within
		for (auto & elem : all_kmers) {
			if (elem.string == minimizer_lookup.front().string)
			{
				found = true;
				break;
			}
		}
		if (!found) { /*if the k-mer is unique, add it*/
			all_kmers.push_back(minimizer_lookup.front());
		}
	}
	/*combining (w,k)-minimizers of a string with (u,k)-end-minimizers for u=1,..,w-1 at both ends of the string,
	all bases of the string will be covered with some minimizer*/
	std::vector<Kmer> end_kmers_left;
	std::vector<Kmer> end_kmers_right;

	for (int u = 1; u < w; ++u) {
		std::vector<Kmer> end_minimizer_lookup_left;
		std::vector<Kmer> end_minimizer_lookup_right;
		/* cycle over the window k+u-1 till k-mers of length, k, can still be made */
		for (int i = 0; i < u; i++) {
			end_minimizer_lookup_left.push_back(Kmer(sequence->genomeString.substr(i, k), i, sequence->identifier));
			end_minimizer_lookup_right.push_back(Kmer(sequence->genomeString.substr(n - k - i, k), n - k - i, sequence->identifier));
		}

		std::sort(end_minimizer_lookup_left.begin(), end_minimizer_lookup_left.end(), sort_Kmers);
		bool found_left = false;
		// Iterate over all minimizers at the left end, add the minimizer of given window only if it is unique
		for (auto & elem : end_kmers_left) {
			if (elem.string == end_minimizer_lookup_left.front().string)
			{
				found_left = true;
				break;
			}
		}
		if (!found_left) {
			end_kmers_left.push_back(end_minimizer_lookup_left.front());
		}

		std::sort(end_minimizer_lookup_right.begin(), end_minimizer_lookup_right.end(), sort_Kmers);
		bool found_right = false;
		// Iterate over all minimizers at the rigth end, add the minimizer of given window only if it is unique
		for (auto & elem : end_kmers_right) {
			if (elem.string == end_minimizer_lookup_right.front().string)
			{
				found_right = true;
				break;
			}
		}
		if (!found_right) {
			end_kmers_right.push_back(end_minimizer_lookup_right.front());
		}
	}

	/*add only unique minimizers from end-minimizers to all_minimizers*/
	
	// Iterate over all minimizers at the left end, add the minimizer to all_kmers only if it is unique
	for (auto & elem1 : end_kmers_left) {
		bool found_left = false;
		for (auto & elem : all_kmers) {
			if (elem.string == elem1.string)
			{
				found_left = true;
				break;
			}
		}
		if (!found_left) {
			all_kmers.push_back(elem1);
		}
	}
	// Iterate over all minimizers at the right end, add the minimizer to all_kmers only if it is unique
	for (auto & elem2 : end_kmers_right) {
		bool found_right = false;
		for (auto & elem : all_kmers) {
			if (elem.string == elem2.string)
			{
				found_right = true;
				break;
			}
		}
		if (!found_right) {
			all_kmers.push_back(elem2);
		}
	}
	std::sort(all_kmers.begin(), all_kmers.end(), sort_pozicija);
	return all_kmers;
};


Kmer_extraction::Kmer_extraction(int x, int y) {
	w = x;
	k = y;
	window_size = w + k - 1;
}
