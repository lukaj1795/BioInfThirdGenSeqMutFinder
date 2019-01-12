#include"Genome.h"
#include"Kmer_extraction.h"
#include"Kmer.h"
#include "FASTAUtility.h"
#include "Alignment.h"
#include "mapping.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <map>
#include "MutationFinder.h"
#include <algorithm>

using namespace std;
using namespace FASTAUtility;
using namespace mapping;


//simple test for mutation finder
void mutation_test(){
	auto g1 = Genome(1);
	auto g2 = Genome(2);

	g1.genomeString = "AACCTTGGAACCG"; //needs to be enlarged so strings have atleast 5 exact kmers for map_to_reference positive match
	g2.genomeString = "AACCTTGGACTCG";

	auto kextraxt = Kmer_extraction(5, 5);
	auto kmer1 = kextraxt.extract(&g1);
	auto kmer2 = kextraxt.extract(&g2);
	int pos = mapping::map_to_reference(&kmer1, &kmer2); //alternative mapping is obsolete, map to reference requires atleast 5 exact kmers within sequences to find position
	std::cout<< pos << "\n";
	std::map<int, std::vector<MutationFinder::MutationOutput>> map;
	for (int i = 0; i < (int)kmer2.size(); i++){
		auto ali = Alignment::Align(kmer1[i].string, kmer2[i].string);
		MutationFinder::map_mutations(kmer1[i].position+ali.first, ali.second, map);

		if (ali.second.size()>0){
			std::cout<< ali.second[0] << "\n";
			std::cout<< ali.second[1] << "\n";
			}
		}

	auto m = MutationFinder::MapToVector(map);
	MutationFinder::output_to_file("mut.csv", m);
}



int main(int argc, char** argv){

	//mutation_test();
	int w = 11;
	int k = 14;
	auto start = chrono::high_resolution_clock::now();
	std::map<int, std::vector<MutationFinder::MutationOutput>> map;
	Alignment::init_vectors(k);
	//reading fasta files
	vector<string> V;
	vector<string> references;	
	std::string name="";
	if (argc >= 3){
		ifstream file(argv[1]);
		ifstream file2(argv[2]);
		name = argv[1];
		if (file.is_open()){
			references = ReadFasta(file);
			std::cout<< "number of references:	" << references.size() << "\n";

			}
		if (file2.is_open()){
			V = ReadFasta(file2);
			std::cout<< "number of sequences:	" << V.size() << "\n";

			}
		file.close();
		file2.close();
		}
	if (argc == 5) {
		w = std::atoi(argv[3]);
		k = std::atoi(argv[4]);
	
	}
	auto finish = chrono::high_resolution_clock::now();
	std::cout<< "Reading input files:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n\n";

	Genome g1 = Genome(0); //reference has identifier 0

	start = chrono::high_resolution_clock::now();
	g1.genomeString.reserve(references[0].length());
	g1.genomeString = std::move(references[0]);

	Kmer_extraction kmer = Kmer_extraction(w, k);
	auto kmer0 = kmer.extract(&g1); //reference k-mer extraction
	finish = chrono::high_resolution_clock::now();
	std::cout<< "k_mer extraction reference:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";
	auto len = g1.genomeString.length();

	g1.genomeString.clear();
	start = chrono::high_resolution_clock::now();
	std::map<int, std::vector<Kmer>> sequence_kmers; //keeping all the k-mers of sequences for mapping, int is identifier
	int counter = 0; //sequence counter for identification
	std::unordered_map<int, int> kmer_position;
	std::unordered_map<int, int> kmer_index;
	std::unordered_map<int, int> pos_index;
	for (int i = 0;i < kmer0.size();i++) {
		kmer_position[kmer0[i].ordering_number_for_string] = kmer0[i].position;
		kmer_index[kmer0[i].ordering_number_for_string] = i;
		pos_index[kmer0[i].position] = i;
	}

	for (auto const& i : V) {
		counter++;
		if (counter % 500 == 0) { //checking progress
			finish = chrono::high_resolution_clock::now();
			std::cout<< "\nstring:	" << counter <<" "<<std::chrono::duration_cast<chrono::seconds>(finish - start).count()  << " s\n";
		}

		Genome g = Genome(counter); //creating new sequence with identifier=counter
		g.genomeString.reserve(i.size());
		g.genomeString = std::move(i);
		auto kmerx = kmer.extract(&g); //extracting k-mers for sequence
		int position = mapping::alternative_mapping(kmer_position, kmerx); //trying to find position of sequence inside the reference
		if (position != -99999) { //if position not equal -99999 position is found
			sequence_kmers.insert(std::pair<int, std::vector<Kmer>>(position, kmerx)); //save position for mutation finding, position isn't 100% accurate
			continue; //don't search the reverse if sequence is mapped
		}
		// REVERSE COMPLEMENT NEEDED, some sequences are in reverse
		Genome gr = Genome(counter);
		auto ir = i;
		reverse(ir.begin(), ir.end()); //creating reverse string
		gr.genomeString.reserve(ir.length());
		gr.genomeString = std::move(ir);
		auto kmery = kmer.extract_complement(&gr); //extracting k-mers of strings complement!!
		int positiony = mapping::alternative_mapping(kmer_position, kmery); //trying to find position of reverse complement
		if (positiony != -99999) { //if position is find
			sequence_kmers.insert(std::pair<int, std::vector<Kmer>>(positiony, kmery)); //save the postion for mutation finding, position isn't 100% accurate
		}
	}

	V.clear();
	int shift = 0; //moving the sequence around the mapped position in the reference
	int flag = 0;
	//int print = 0;
	for (auto seq_K : sequence_kmers) { //for every mapped sequence

		auto pos = seq_K.first; //save the position to which it is aligned
		if (pos < 0) {
			continue;
		}
		auto kmers = seq_K.second; //k-mers of mapped sequence
		shift = 0; //reset the shift for current sequence
		flag = 0;
		//std::cout<< "rad na kmerima sekv " << print << "\n";
		//print++;
		for (int i = pos_index[pos]; i < kmer0.size(); i++) {
			if (kmer0[i].position + 15 > pos) { //start 15 positions earlier, mapping is not 100% accurate
				flag++;
				if (flag > kmers.back().position + kmers.size()/20) //reference has been searched enough for this current sequence
					break;
				auto kmer_ref = kmer0[i];
				for (int j = shift; j < shift + kmers.size()/10; j++) {
					if (len > kmer0[i].position + kmers[j].position) { //if we are inside the reference
						if (mapping::check_match(kmer_ref.string, kmers[j].string)) { //checking if referent k-mer and sequence k-mer are similar enough to enter the alignment process
							auto aligned = Alignment::Align(kmer_ref.string, kmers[j].string); //alignment of two k-mers
							if (!aligned.second.empty()) { //if only 1 mutation found, possible outcome for the each base : MATCH; INSERT; DELETE; SUBSTITUTION
								MutationFinder::map_mutations(kmer_ref.position + aligned.first, aligned.second, map); //keep the alignment data
							}
						}
					}
				}
				if (flag > kmers.size()/20) {
					shift++; //when some amount of k-mers are checked for mutations, start shifting the sequence along the reference
					if (shift + kmers.size() / 10 == kmers.size()) { //when the end of sequence is reached stop shifting
						shift--; //shift stop
					}
				}
			}
		}
	}
	sequence_kmers.clear();
	kmer0.clear();

	auto m = MutationFinder::MapToVector(map); //extraction of most often mutation (or MATCH) for each position of reference

	map.clear();

	//print - csv
	MutationFinder::output_to_file(std::to_string(w)+"w_"+std::to_string(k)+"k_"+name.substr(0,name.length()-6)+"_mut.csv", m);

	finish = chrono::high_resolution_clock::now();
	std::cout<< "Total:	" << std::chrono::duration_cast<chrono::milliseconds>(finish - start).count()  << " ms\n";

	return 0;
}

