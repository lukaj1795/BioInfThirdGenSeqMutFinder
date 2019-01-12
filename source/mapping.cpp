#include<string>
#include "Kmer.h"
#include "Alignment.h"
#include "mapping.h"
#include <iostream>
#include <map>
#include <set>
#include "Kmer_mapping.h"
#include <unordered_map>


namespace mapping {

	namespace {
		const int error = -99999;
		const int mutation_number = 2;
		/*std::map<char, int> k1_count = {
			{ 'C', 0 },
			{ 'A', 0 },
			{ 'T', 0 },
			{ 'G', 0 }
		};
		std::map<char, int> k2_count = {
			{ 'C', 0 },
			{ 'A', 0 },
			{ 'T', 0 },
			{ 'G', 0 }
		};*/
	}
	//std::vector<char> bases = { 'A', 'C', 'T', 'G' };

	int map_to_reference(std::vector<Kmer>* reference, std::vector<Kmer>* sequence)
	{

		int kmer_size = reference->front().string.size();
		int n = 12;
		int match = 0;
		int flag1 = 0;
		int counter = 0;
		std::vector<int> seq_pos_match;
		std::vector<int> ref_pos;
		std::vector<int> seq_pos;
		std::map <int,Kmer_mapping> matches;
		//std::vector<Kmer>::iterator speed_up;
		counter = 0;
		int step = 30;
		bool mapped;
	
		for (auto kmer_ref_iterator = reference->begin(); kmer_ref_iterator != reference->end(); ++kmer_ref_iterator) {
		
			counter++;
			if (match == 0) { //skip 30 kmers until atleast 1 absolute match is found
				if (counter < step) {
					continue;
				}
				counter = 0; //reset counter for skipping STEP number of references
			}
			if (std::next(kmer_ref_iterator, 1) == reference->end() || std::next(kmer_ref_iterator, 2) == reference->end() || std::next(kmer_ref_iterator, 3) == reference->end() || std::next(kmer_ref_iterator, 4) == reference->end() || std::next(kmer_ref_iterator, 5) == reference->end() || std::next(kmer_ref_iterator, 6) == reference->end() || std::next(kmer_ref_iterator, 7) == reference->end() || std::next(kmer_ref_iterator, 8) == reference->end() || std::next(kmer_ref_iterator, 9) == reference->end() || std::next(kmer_ref_iterator, 10) == reference->end()) {
				return -99999; //if the end of reference is reached
			}
			auto kmer_ref = *kmer_ref_iterator; //check few values (with some step)
			auto kmer_ref2_it = std::next(kmer_ref_iterator, 2);
			auto kmer_ref4_it = std::next(kmer_ref_iterator, 4);
			auto kmer_ref6_it = std::next(kmer_ref_iterator, 6);

			Kmer kmer_ref2 = *kmer_ref2_it;
			Kmer kmer_ref4 = *kmer_ref6_it;
			Kmer kmer_ref6 = *kmer_ref6_it;

			flag1 = 0;
			//iterate through sequence
			for (auto kmer_seq_iterator = sequence->begin(); kmer_seq_iterator != sequence->end(); ++kmer_seq_iterator) {
				auto kmer_seq = *kmer_seq_iterator;
				flag1++;
				if (flag1 < 3) { //skip first few k-mers, since they are probably end-mers and different from sequence k-mers
					continue;
				}
				if (flag1 > 209) { //if match isn't found in the first 209 k-mers, stop search
					break;
				}
				if (match == 0) //untill 1 match is found, check few k-mers at once with big step between searches
					mapped = (kmer_ref.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref2.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref4.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref6.ordering_number_for_string == kmer_seq.ordering_number_for_string); //|| kmer_ref8.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref10.ordering_number_for_string == kmer_seq.ordering_number_for_string);
				
				else //after 1 match is found, search one by one
					mapped = kmer_ref.ordering_number_for_string == kmer_seq.ordering_number_for_string;
				if(mapped){
					counter++;
					match++; //match count
					if (match > 7) {
						return kmer_ref.position - kmer_seq.position; //if enough matches found, return mapped position
						//break;
					}
					//OBSOLETE CODE, reason: TOO SLOW, but better precision
					/*if(matches.find(kmer_seq.position) == matches.end())
					{// key kmer_seq.position doesn't exist
						Kmer_mapping new_position = Kmer_mapping(kmer_seq.position, kmer_ref.position);
						auto result = matches.insert(std::pair<int, Kmer_mapping>(kmer_seq.position, new_position));
						seq_pos_match.push_back(kmer_seq.position);
					}
					else { //key already exists.. add new reference position to Kmer_mapping
						matches.at(kmer_seq.position).add_position(kmer_ref.position);
					}*/

				}
			}
		}
		/*if (match <= 5) { //if sequence is very damaged
			return -99999;
		}
		//check if the 3 sequal matches exist, map only quality sequnces
		for (auto it1 = seq_pos_match.begin(); it1 != seq_pos_match.end(); ++it1) {
			if (std::next(it1, 1) != seq_pos_match.end()) {
				auto it2 = std::next(it1, 1);
				if (std::next(it1, 2) != seq_pos_match.end()) {
					auto it3 = std::next(it1, 2);
					if (std::next(it1, 3) != seq_pos_match.end()) {
						auto a = matches.at(*it1).reference_positions;
						auto b = matches.at(*it2).reference_positions;
						auto c = matches.at(*it3).reference_positions;
						for (auto aa : a) {
							for (auto bb : b) {
								if (bb - aa >= *it2 - *it1 - n && bb - aa <= *it2 - *it1 + n) {
									for (auto cc : c) {
										if (cc - bb >= *it3 - *it2 - n && cc - bb <= *it3 - *it2 + n) {
											return(aa - *it1);
										}
									}
								}
							}
						}
					}
				}
			}
		}*/
		//std::cout << "\nNema PREKLAPANJa!!	" << match;
		return -99999;
	}

	int alternative_mapping(const std::unordered_map<int,int> &reference, const std::vector<Kmer> & sequence) {

		int match = 0;

		int counter = 0;
		int threshold = sequence[0].string.length();


		for (auto seq = sequence.begin(); seq != sequence.end(); ++seq) {
			auto seq_K = *seq;
			auto ss = seq_K.ordering_number_for_string;
			if (reference.find(ss)!=reference.end()) {
				///std::cout << ref.string << " " << seq_K.string << "\n";
				return reference.at(ss) - seq_K.position;
			}
			if (counter > threshold)
			{
				break;
			}
			counter++;
		}

		return error;
	}
	
	bool check_match(const std::string &k1, const std::string &k2) {

		///maybe unnecessary
		/*if (k1.ordering_number_for_string == k2.ordering_number_for_string){
			return true;
			}*/
		/*for (auto a : k1_count) {
			
			k1_count.at(a.first) = 0;
			//k2_count.at(a.first) = 0;
		}

		//int j = 0;

		for (std::string::size_type i = 0; i < k1.string.size(); ++i) {
			k1_count.at(k1.string[i])++;
			k1_count.at(k2.string[i])--;
		}
		//return diff < mutation_number * 2;
		//finished counting;
		return !((abs(k1_count['A']) + abs(k1_count['T']) + abs(k1_count['C']) + abs(k1_count['G'])) > mutation_number * 2);
		*/
		int diff = 0;
		for (int i = 0;i < k1.size();i++) {
			diff += k1[i] != k2[i];
		}
		return diff < mutation_number * 2;
		
		}

	}
	
	
