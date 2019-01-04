#include"Kmer_extraction.h"
#include"Genome.h"
#include"Kmer.h"
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>

bool sort_Kmers(Kmer k1, Kmer k2)
{
	return (k1.ordering_number_for_string < k2.ordering_number_for_string);
}
bool sort_pozicija(Kmer k1, Kmer k2)
{
	return (k1.position < k2.position);
}
std::vector<Kmer> Kmer_extraction::extract(Genome *sequence) {

	int n = sequence->genomeString.size();
	std::map<int,Kmer> all_kmers;
	std::set<Kmer> minimizer_lookup;
	std::vector<Kmer> all_return;
	/*the length of String*/
	for (int j = 0; j <= n - window_size; ++j) {
		/* cycle over the window till k-mers of length, k, can still be made */
		for (int i = j; i < j + w; i = i + 1) {
			minimizer_lookup.insert(Kmer(std::move(sequence->genomeString.substr(i, k)), i, sequence->identifier));
		}
		/*only smallest k-mer will become minimizer of the given window*/
		auto someElementIterator = minimizer_lookup.begin();
		auto value = *someElementIterator;
		all_kmers.insert(std::pair<int, Kmer>(value.position,value));
		minimizer_lookup.clear();
	}
	/*combining (w,k)-minimizers of a string with (u,k)-end-minimizers for u=1,..,w-1 at both ends of the string,
	all bases of the string will be covered with some minimizer*/
	std::set<Kmer> end_kmers_left;
	std::set<Kmer> end_kmers_right;

	std::set<Kmer> end_minimizer_lookup_left;
	std::set<Kmer> end_minimizer_lookup_right;

	for (int u = 1; u < w; ++u) {
		
		/* cycle over the window k+u-1 till k-mers of length, k, can still be made */
		for (int i = 0; i < u; i++) {
			end_minimizer_lookup_left.insert(Kmer(sequence->genomeString.substr(i, k), i, sequence->identifier));
			end_minimizer_lookup_right.insert(Kmer(sequence->genomeString.substr(n - k - i, k), n - k - i, sequence->identifier));
		}

		// Iterate over all minimizers at the left end, add the minimizer of given window only if it is unique
		auto someElementIterator =end_minimizer_lookup_left.begin();
		auto value = *someElementIterator;
		end_kmers_left.insert(value);
		someElementIterator = end_minimizer_lookup_right.begin();
		value = *someElementIterator;
		end_kmers_right.insert(value);
		end_minimizer_lookup_left.clear();
		end_minimizer_lookup_right.clear();
	}

	/*add only unique minimizers from end-minimizers to all_minimizers*/
	
	for (auto & elem : end_kmers_left) {
		all_kmers.insert(std::pair<int, Kmer>(elem.position, elem));
	}
	for (auto & elem : end_kmers_right) {
		all_kmers.insert(std::pair<int, Kmer>(elem.position, elem));
	}
	all_return.reserve(all_kmers.size());
	for (auto i : all_kmers) {
		all_return.push_back(std::move(i.second));
	}
	return all_return;
};


Kmer_extraction::Kmer_extraction(int x, int y) {
	w = x;
	k = y;
	window_size = w + k - 1;
}

std::vector<Kmer> Kmer_extraction::extract_complement(Genome *sequence) {
	std::string complement= "";
	std::map<char, char> komplement = {
		{ 'C', 'G' },
		{ 'A', 'T' },
		{ 'T', 'A' },
		{ 'G', 'C' }
	};
	for (auto i : sequence->genomeString) {
		complement += komplement.at(i);
	}
	Genome g = Genome(sequence->identifier);
	g.genomeString = complement;
	return Kmer_extraction::extract(&g);
}