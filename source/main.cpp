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

using namespace std;
using namespace FASTAUtility;

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
	std::cout << pos << "\n";
	std::map<int, std::vector<MutationFinder::MutationOutput>> map;
	for (int i = 0; i < (int)kmer2.size(); i++){
		auto ali = Alignment::Align(kmer1[i], kmer2[i]);
		MutationFinder::map_mutations(kmer1[i].position+ali.first, ali.second, map);

		if (ali.second.size()>0){
			std::cout << ali.second[0] << "\n";
			std::cout << ali.second[1] << "\n";
			}
		}

	auto m = MutationFinder::MapToVector(map);
	MutationFinder::output_to_file("mut.csv", m);

	int x;
	cin >> x;
	}



void main(int argc, char** argv){

	//mutation_test();
	int w = 20;
	int k = w;
	auto start = chrono::high_resolution_clock::now();
	std::map<int, std::vector<MutationFinder::MutationOutput>> map;
	//reading fasta files
	vector<string> V;
	vector<string> references;	
	std::string name="";
	if (argc == 3){
		ifstream file(argv[1]);
		ifstream file2(argv[2]);
		name = argv[1];
		if (file.is_open()){
			references = ReadFasta(file);
			cout << "number of references:	" << references.size() << "\n";

			}
		if (file2.is_open()){
			V = ReadFasta(file2);
			cout << "number of sequences:	" << V.size() << "\n";

			}
		file.close();
		file2.close();
		}
	auto finish = chrono::high_resolution_clock::now();
	cout << "Reading input files:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n\n";

	Genome g1 = Genome(0); //reference has identifier 0
	//g1.genomeString = "GTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGATGCATCAGCGAGCAGCGATCGTAAGACCGCTAGCTAAGCTCGGCCTACGATAACGACATCAGCTACGATGCATCGATCTGATCGAGCATGCTGAGCAGCGTACTATGCGTAGTCATGCTGAGTGTCTTGGTCAGCAAAATGCATCGATCGACATGGTGTTCGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGAGGTCGTCACGTTACTCACAAGCAT";
	//g2.genomeString = "CGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGA";
	start = chrono::high_resolution_clock::now();
	g1.genomeString = references[0];

	Kmer_extraction kmer = Kmer_extraction(w, k);
	auto kmer0 = kmer.extract(&g1);
	finish = chrono::high_resolution_clock::now();
	cout << "k_mer extraction reference:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";

	start = chrono::high_resolution_clock::now();
	std::map<int, int> pokrivenost1;
	std::map<int, std::vector<Kmer>> sequence_kmers; //keeping all the kmers, int is identifier
	int counter = 0;
	for (auto const& i : V) {
		counter++;
		//cout << "\nstring:	" << counter;
		Genome g = Genome(counter);
		g.genomeString = i;
		//auto start1 = chrono::high_resolution_clock::now();
		auto kmerx = kmer.extract(&g);
		//auto finish1 = chrono::high_resolution_clock::now();
		//cout << "k_mer extraction seq:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish1 - start1).count() / 1e6 << " ms\n";
		//start1 = chrono::high_resolution_clock::now();
		int position = mapping::map_to_reference(&kmer0, &kmerx);
		//finish1 = chrono::high_resolution_clock::now();
		//cout << "mapping:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish1 - start1).count() / 1e6 << " ms\n";
		if (position != -99999) {
			cout <<"string:	"<< counter << "\nNORMAL preklapanje:	" << position << "	+-20 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int, int>(position, counter));
			sequence_kmers.insert(std::pair<int, std::vector<Kmer>>(position, kmerx));
			continue; //don't search the reverse if it is mapped
			///get aligned strings
			///map mutations
		}
		// REVERSE NEEDED, some sequences are in reverse
		Genome gr = Genome(counter);
		auto ir = i;
		reverse(ir.begin(), ir.end());
		gr.genomeString = ir;
		auto kmery = kmer.extract_complement(&gr);
		int positiony = mapping::map_to_reference(&kmer0, &kmery);
		if (positiony != -99999) {
			cout <<"string:	"<< counter << "\nREVERSE  preklapanje:	" << positiony << "	+-20 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int, int>(positiony, counter));
			continue;
			///get aligned strings
			///map mutations
		}
		/* COMPLEMENT FOUND NOT IMPORTANT AGAIN
		auto kmerz = kmer.extract_complement(&gr);
		int positionz = mapping::map_to_reference(&kmer0, &kmerz);
		if (positionz != -99999) {
			cout << "string:	" << counter << "\nREVERSE complement preklapanje:	" << positionz << "	+-20 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int, int>(positionz, counter));
			///get aligned strings
			///map mutations
		}
		auto kmerw = kmer.extract_complement(&g);
		int positionw = mapping::map_to_reference(&kmer0, &kmerw);
		if (positiony != -99999) {
			cout << "string:	" << counter << "\nNORMAL complement preklapanje:	" << positionw << "	+-20 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int, int>(positionw, counter));
			///get aligned strings
			///map mutations
		}*/
	}


	for (auto seq_K : sequence_kmers){
		
		auto pos = seq_K.first;
		auto kmers = seq_K.second;


		for (int i = 0; i < kmer0.size(); i++){
			if (kmer0[i].position > pos){
				for (int j = 0; j < kmers.size();j++){
					auto kmer_ref = kmer0[i - 1 + j];
					///std::cerr << "\n";
					///std::cerr << kmer_ref.string << "\n";
					///std::cerr << kmers[j].string << "\n";
					auto aligned = Alignment::Align(kmer_ref,kmers[j]);
					/*if (aligned.second.size()){
						std::cerr << aligned.second[0] << "\n";
						std::cerr << aligned.second[1] << "\n";
						}*/
					MutationFinder::map_mutations(kmer_ref.position+aligned.first,aligned.second,map);
					}
				break;
				}
			
			}

		
		}



	auto m = MutationFinder::MapToVector(map);

	MutationFinder::output_to_file(name.substr(0,name.length()-6)+"_mut.csv", m);


	finish = chrono::high_resolution_clock::now();
	cout << "Alignment:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";


	for (auto i : pokrivenost1) {			
		cout << "\nPosition:	" << i.first << "	string:	" << i.second;
	}

		int x;
		cin >> x;

		///mutation_test();

}

