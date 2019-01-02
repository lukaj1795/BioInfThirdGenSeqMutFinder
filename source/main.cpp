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

	g1.genomeString = "AACCTTGGAACCG";
	g2.genomeString = "AACCTTGGACTCG";

	auto kextraxt = Kmer_extraction(5, 5);
	auto kmer1 = kextraxt.extract(&g1);
	auto kmer2 = kextraxt.extract(&g2);
	int pos = mapping::alternative_mapping(&kmer1, &kmer2);
	std::cout << pos << "\n";
	std::map<int, std::vector<MutationFinder::MutationOutput>> map;
	for (int i = 0; i < kmer2.size(); i++){
		auto ali = Alignment::Align(kmer1[i], kmer2[i]);
		MutationFinder::map_mutations(kmer1[i].position, ali, map);

		if (ali.size()>0){
			std::cout << ali[0] << "\n";
			std::cout << ali[1] << "\n";
			}
		}

	auto m = MutationFinder::MapToVector(map);
	MutationFinder::output_to_file("mut.csv", m);

	int x;
	cin >> x;
	}



void main(int argc, char** argv){

	mutation_test();

	int k = 20;
	auto start = chrono::high_resolution_clock::now();
	std::map<int, std::vector<MutationFinder::MutationOutput>> map;
	//reading fasta files
	vector<string> V;
	vector<string> references;
	if (argc == 3){
		ifstream file(argv[1]);
		ifstream file2(argv[2]);
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
	//cout << reference << "\n";
	//cout << V[0] << "\n";


	Genome g1 = Genome(1);
	Genome g2 = Genome(2);
	//g1.genomeString = "GTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGATGCATCAGCGAGCAGCGATCGTAAGACCGCTAGCTAAGCTCGGCCTACGATAACGACATCAGCTACGATGCATCGATCTGATCGAGCATGCTGAGCAGCGTACTATGCGTAGTCATGCTGAGTGTCTTGGTCAGCAAAATGCATCGATCGACATGGTGTTCGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGAGGTCGTCACGTTACTCACAAGCAT";
	//g2.genomeString = "CGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGA";
	start = chrono::high_resolution_clock::now();
	g1.genomeString = references[0];
	//g2.genomeString = V[0];

	Kmer_extraction kmer = Kmer_extraction(20, k);
	auto kmer1 = kmer.extract(&g1);
	finish = chrono::high_resolution_clock::now();
	cout << "k_mer extraction reference:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";

	start = chrono::high_resolution_clock::now();
	std::map<int, int> pokrivenost1;
	int counter = 0;
	for (auto i : V) {
		counter++;
		///cout << "\n\nstring:	" << counter << "	duljina:	" << i.size();
		g2.genomeString = i;
		auto kmer2 = kmer.extract(&g2);
		int position = mapping::map_to_reference(&kmer1, &kmer2);
		if (position != -99999) {
			//cout << "\npozicija unutar reference gdje bi se trebali poceti preklapati:	" << position << "	+-10 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int, int>(position, counter));
			///get aligned strings
			///map mutations

			}




		auto m = MutationFinder::MapToVector(map);
		MutationFinder::output_to_file("lambda_mut.csv", m);


		finish = chrono::high_resolution_clock::now();
		cout << "Alignment:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";


		for (auto i : pokrivenost1) {
			cout << "\nPosition:	" << i.first << "	string:	" << i.second;

			}



		/*
		cout << "\n\n" << "string 1: " << g1.genomeString << "\n\n";
		for (auto i : kmer1)
		std::cout <<"kmer: "<< i.string << " vrijednost: " << i.ordering_number_for_string<<" pozicija: " << i.position << " sekvenca: " << i.identifier <<"\n";

		cout <<"\n\n"<< "string 2: " << g2.genomeString << "\n\n";
		for (auto i : kmer2)
		std::cout << "kmer: " << i.string << " vrijednost: " << i.ordering_number_for_string << " pozicija: " << i.position << " sekvenca: " << i.identifier << "\n";
		*/


		//int x;
		//cin >> x;



		///mutation_test();

		}

	}