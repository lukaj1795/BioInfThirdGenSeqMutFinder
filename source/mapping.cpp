#include<string>
#include "Kmer.h"
#include "Alignment.h"
#include "mapping.h"
#include <iostream>
#include <map>
#include <set>
#include "Kmer_mapping.h"

namespace mapping {
	int error = -99999;
	int map_to_reference(std::vector<Kmer>* reference, std::vector<Kmer>* sequence)
	{
		int n = 10;
		int match = 0;
		int flag1 = 0;
		int counter = 0;
		std::vector<int> seq_pos_match;
		std::vector<int> ref_pos;
		std::vector<int> seq_pos;
		std::map <int,Kmer_mapping> matches;
		std::vector<Kmer>::iterator speed_up;
		
		for (auto kmer_ref_iterator = reference->begin(); kmer_ref_iterator != reference->end(); ++kmer_ref_iterator) {
			
			auto kmer_ref = *kmer_ref_iterator;

			if (match>12) { //positive id on match--we can now find starting point
				break;
			}
			counter = 0;
			flag1 = 0;
			for (auto kmer_seq_iterator = sequence->begin(); kmer_seq_iterator != sequence->end(); ++kmer_seq_iterator) {
				auto kmer_seq = *kmer_seq_iterator;
				flag1++;
				if (flag1 > 39) {
					break;
				}
				if (!mapping::check_match(kmer_ref, kmer_seq)) {
					continue;
				}
				auto mapped = Alignment::Align(kmer_ref, kmer_seq);
				if (!mapped.empty()) { //VERY IMPORTANT!!! ALL THE VALUES OF ALREADY ALIGNED KMERS GET BACK HERE INSIDE VECTOR MAPPED- CAN BE USED FOR MUTATION FINDING
					//std::cout << "\nPOCETAKref:	" << kmer_ref.string << "\nPrviKmerSeq:	" << kmer_seq.string << "\n";
					match++;
					counter++;
					
					if(matches.find(kmer_seq.position) == matches.end())
					{// key kmer_seq.position doesn't exist
						Kmer_mapping new_position = Kmer_mapping(kmer_seq.position, kmer_ref.position);
						auto result = matches.insert(std::pair<int, Kmer_mapping>(kmer_seq.position, new_position));
						seq_pos_match.push_back(kmer_seq.position);
					}
					else { //key already exists.. add new reference position to Kmer_mapping
						matches.at(kmer_seq.position).add_position(kmer_ref.position);
					}
					if (counter > 2) {
						break;
					}
				}
			}
		}
		//for (auto &i : seq_pos_match) { //al the positions in sequence file where string appears
		for (auto it1 = seq_pos_match.begin(); it1 != seq_pos_match.end(); ++it1) {
			if (std::next(it1, 1) != seq_pos_match.end()) {
				auto it2 = std::next(it1, 1);
				if (std::next(it1, 2) != seq_pos_match.end()) {
					auto it3 = std::next(it1, 2);
					if (std::next(it1, 3) != seq_pos_match.end()) {
						auto it4 = std::next(it1, 3);
						if (std::next(it1, 4) != seq_pos_match.end()) {
							auto it5 = std::next(it1, 4);
							if (std::next(it1, 5) != seq_pos_match.end()) {
								auto it6 = std::next(it1, 5);

								auto a = matches.at(*it1).reference_positions;
								auto b = matches.at(*it2).reference_positions;
								auto c = matches.at(*it3).reference_positions;
								auto d = matches.at(*it4).reference_positions;
								auto e = matches.at(*it5).reference_positions;
								//auto f = matches.at(*it6).reference_positions;

								for (auto aa : a) {
									
									for (auto bb : b) {
										if (bb - aa >= *it2 - *it1 - n && bb - aa <= *it2 - *it1 + n) {
											for (auto cc : c) {
												if (cc - bb >= *it3 - *it2 - n && cc - bb <= *it3 - *it2 + n) {
													for (auto dd : d) {
														if (dd - cc >= *it4 - *it3 - n && dd - cc <= *it4 - *it3 + n) {
															for (auto ee : e) {
																if (ee - dd >= *it5 - *it4 - n && ee - dd <= *it5 - *it4 + n) {
																	//for (auto ff : f) {
																		//if (ff - ee >= *it6 - *it5 - n && ff - ee <= *it6 - *it5 + n) {
																			std::cout << "\nNADENO PREKLAPANJE!!	"<< match;
																			return(aa - *it1);
																			
																		//}
																	//}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//std::cout << "\nNema PREKLAPANJa!!	" << match;
		return -99999;
	}
	int alternative_mapping(std::vector<Kmer>* reference, std::vector<Kmer>* sequence) {

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
	}
	bool mapping::check_match(Kmer k1, Kmer k2) {
		std::map<char, int> k1_count = {
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
		};
		int j = 0;
		for (std::string::size_type i = 0; i < k1.string.size(); ++i) {
			j=k1_count.at(k1.string[i]);
			k1_count.at(k1.string[i]) = ++j;
			j=k2_count.at(k2.string[i]);
			k2_count.at(k2.string[i]) = ++j;
		}
		//finished counting;
		if (abs(k1_count.at('A') - k2_count.at('A')) + abs(k1_count.at('T') - k2_count.at('T')) + abs(k1_count.at('C') - k2_count.at('C')) + abs(k1_count.at('G') - k2_count.at('G')) > 3) {
			return false;
		}
		return true;
	}
}