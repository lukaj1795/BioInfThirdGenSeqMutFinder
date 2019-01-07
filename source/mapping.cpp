#include<string>
#include "Kmer.h"
#include "Alignment.h"
#include "mapping.h"
#include <iostream>
#include <map>
#include <set>
#include "Kmer_mapping.h"


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
		int step = 100;
		bool mapped;
	
		for (auto kmer_ref_iterator = reference->begin(); kmer_ref_iterator != reference->end(); ++kmer_ref_iterator) {
		
			counter++;
			if (match == 0) {
				if (counter < step) {
					continue;
				}
				counter = 0; //reset counter for skipping STEP number of references
			}
			if (std::next(kmer_ref_iterator, 1) == reference->end() || std::next(kmer_ref_iterator, 2) == reference->end() || std::next(kmer_ref_iterator, 3) == reference->end() || std::next(kmer_ref_iterator, 4) == reference->end() || std::next(kmer_ref_iterator, 5) == reference->end() || std::next(kmer_ref_iterator, 6) == reference->end() || std::next(kmer_ref_iterator, 7) == reference->end() || std::next(kmer_ref_iterator, 8) == reference->end() || std::next(kmer_ref_iterator, 9) == reference->end() || std::next(kmer_ref_iterator, 10) == reference->end()) {
				return -99999;
			}
			auto kmer_ref = *kmer_ref_iterator;
			auto kmer_ref2_it = std::next(kmer_ref_iterator, 2);
			auto kmer_ref4_it = std::next(kmer_ref_iterator, 4);
			auto kmer_ref6_it = std::next(kmer_ref_iterator, 6);
			auto kmer_ref8_it = std::next(kmer_ref_iterator, 8);
			auto kmer_ref10_it = std::next(kmer_ref_iterator, 10);
			Kmer kmer_ref2 = *kmer_ref2_it;
			Kmer kmer_ref4 = *kmer_ref6_it;
			Kmer kmer_ref6 = *kmer_ref6_it;
			Kmer kmer_ref8 = *kmer_ref8_it;
			Kmer kmer_ref10 = *kmer_ref10_it;
			flag1 = 0;
			for (auto kmer_seq_iterator = sequence->begin(); kmer_seq_iterator != sequence->end(); ++kmer_seq_iterator) {
				auto kmer_seq = *kmer_seq_iterator;
				flag1++;
				if (flag1 < 3) {
					continue;
				}
				if (flag1 > 209) {
					break;
				}
				//if (!mapping::check_match(kmer_ref, kmer_seq)) {
					//continue;
				//}
				//auto mapped = Alignment::Align_int(kmer_ref, kmer_seq); //HERE IS BIGGEST TIME CONSUMER NEEDS TO BE CHANGED
				if (match == 0)
					mapped = (kmer_ref.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref2.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref4.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref6.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref8.ordering_number_for_string == kmer_seq.ordering_number_for_string || kmer_ref10.ordering_number_for_string == kmer_seq.ordering_number_for_string);
				//if (mapped > (kmer_size * 4 - 4 - 2 * 2)) { //VERY IMPORTANT!!! ALL THE VALUES OF ALREADY ALIGNED KMERS GET BACK HERE INSIDE VECTOR MAPPED- CAN BE USED FOR MUTATION FINDING
				else
					mapped = kmer_ref.ordering_number_for_string == kmer_seq.ordering_number_for_string;
				if(mapped){
					counter++;

					match++;
					if (match > 7) {
						return kmer_ref.position - kmer_seq.position;
						//break;
					}
				
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
		/*if (match <= 5) {
			return -99999;
		}
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

	/*int alternative_mapping(std::vector<Kmer>* reference, std::vector<Kmer>* sequence) {

		int match = 0;

		int counter = 0;
		int threshold = 12;

		for (auto kmer = reference->begin(); kmer != reference->end(); ++kmer) {
			auto ref = *kmer;
			counter = 0;
			for (auto seq = sequence->begin(); seq != sequence->end(); ++seq) {
				auto seq_K = *seq;
				auto sr = ref.ordering_number_for_string;
				auto ss = seq_K.ordering_number_for_string;
				if (sr == ss) {
					///std::cout << ref.string << " " << seq_K.string << "\n";
					return ref.position - seq_K.position;
				}
				if (counter > threshold)
				{
					break;
				}
				counter++;
			}

		}
		return error;
	}*/
	
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
	
	
