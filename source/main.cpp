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

	//mutation_test();

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

	Genome g1 = Genome(0); //reference has identifier 0
	//g1.genomeString = "GTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGATGCATCAGCGAGCAGCGATCGTAAGACCGCTAGCTAAGCTCGGCCTACGATAACGACATCAGCTACGATGCATCGATCTGATCGAGCATGCTGAGCAGCGTACTATGCGTAGTCATGCTGAGTGTCTTGGTCAGCAAAATGCATCGATCGACATGGTGTTCGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGAGGTCGTCACGTTACTCACAAGCAT";
	//g2.genomeString = "CGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGA";
	start = chrono::high_resolution_clock::now();
	g1.genomeString = references[0];

	Kmer_extraction kmer = Kmer_extraction(20, k);
	auto kmer0 = kmer.extract(&g1);
	finish = chrono::high_resolution_clock::now();
	cout << "k_mer extraction reference:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";

	start = chrono::high_resolution_clock::now();
	std::map<int, int> pokrivenost1;
	std::map<int, std::vector<Kmer>> sequence_kmers; //keeping all the kmers, int is identifier
	int counter = 0;
	for (auto i : V) {
		counter++;
		cout << "\nstring:	" << counter;
		Genome g = Genome(counter);
		g.genomeString = i;
		auto kmerx = kmer.extract(&g);
		sequence_kmers.insert(std::pair<int, std::vector<Kmer>>(counter, kmerx));
		int positionx = mapping::map_to_reference(&kmer0, &kmerx);
		if (positionx != -99999) {
			cout << "\npozicija unutar reference gdje bi se trebali poceti preklapati:	" << positionx << "	+-10 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int, int>(positionx, counter));
			///get aligned strings
			///map mutations

		}
		cout << "\ncomplement:	";
		auto kmery = kmer.extract_complement(&g);
		int positiony = mapping::map_to_reference(&kmer0, &kmery);
		if (positiony != -99999) {
			cout << "\npozicija unutar reference gdje bi se trebali poceti preklapati:	" << positiony << "	+-10 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int, int>(positiony, counter));
			///get aligned strings
			///map mutations

		}
	}
		auto m = MutationFinder::MapToVector(map);
		MutationFinder::output_to_file("lambda_mut.csv", m);


		finish = chrono::high_resolution_clock::now();
		cout << "Alignment:	" << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";


		for (auto i : pokrivenost1) {
			cout << "\nPosition:	" << i.first << "	string:	" << i.second;

			}

		int x;
		cin >> x;

		///mutation_test();

}

