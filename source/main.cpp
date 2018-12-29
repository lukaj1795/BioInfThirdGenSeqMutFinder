#include"Genome.h"
#include"Kmer_extraction.h"
#include"Kmer.h"
#include "FASTAUtility.h"
#include "Alignment.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>

using namespace std;
using namespace FASTAUtility;

void main(int argc,char** argv){
	auto start = chrono::high_resolution_clock::now();
	//reading fasta files
	vector<string> V;
	vector<string> references;
	if (argc == 3){
		ifstream file(argv[1]);
		ifstream file2(argv[2]);
		if (file.is_open()){
			references = ReadFasta(file);
			cout << references.size() << "\n";

		}
		if (file2.is_open()){
			V = ReadFasta(file2);
			cout << V.size() << "\n";

			}
		file.close();
		file2.close();
	}
	auto reference = references[0];
	auto finish = chrono::high_resolution_clock::now();
	//cout << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";
	//cout << reference << "\n";
	//cout << V[0] << "\n";
	Genome g1= Genome(1);
	Genome g2 = Genome(2);
///	Genome g3 = Genome(3);
///	Genome g4 = Genome(4);
///	g3.genomeString = "AGGCTGTCACGTCAGTGTGCCAACGTGCAACGCCTGACCTGACTGGGTCACTGACTG";
///	g2.genomeString = "CGTCAGTGTGCCAACGTGCAACGCCT";
	//g1.genomeString = "AACGCCTGACCTGACTGGGTCA";
	//g4.genomeString = "AACGCCTAGCCTGACTGGGTCA";
	start = chrono::high_resolution_clock::now();
	g1.genomeString = V[0];
	g2.genomeString = reference;
	Kmer_extraction kmer = Kmer_extraction(10,10);
	//std::vector<Kmer> svi_kmerovi_sekvence;
	auto kmer1=kmer.extract(&g1);
	auto kmer4=kmer.extract(&g2);

	
	//cout << "string 1: " << g1.genomeString << "\n\n";

	//auto svi_kmerovi_sekvence = kmer.extract(&g1);
	finish = chrono::high_resolution_clock::now();
	
	cout << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";
	start = chrono::high_resolution_clock::now();
	Alignment::Align(kmer1,kmer4);
	finish = chrono::high_resolution_clock::now();
	cout << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";
	//for (auto i : svi_kmerovi_sekvence)
		//std::cout <<"kmer: "<< i.string << " vrijednost: " << i.ordering_number_for_string<<" pozicija: " << i.position << " sekvenca: " << i.identifier <<"\n";

	//cout <<"\n\n"<< "string 2: " << g2.genomeString << "\n\n";
	/*svi_kmerovi_sekvence = kmer.extract(&g2);
	for (auto i : svi_kmerovi_sekvence)
		std::cout << "kmer: " << i.string << " vrijednost: " << i.ordering_number_for_string << " pozicija: " << i.position << " sekvenca: " << i.identifier << "\n";

	cout <<"\n\n"<< "string 3: " << g3.genomeString << "\n\n";
	svi_kmerovi_sekvence = kmer.extract(&g3);
	for (auto i : svi_kmerovi_sekvence)
		std::cout << "kmer: " << i.string << " vrijednost: " << i.ordering_number_for_string << " pozicija: " << i.position << " sekvenca: " << i.identifier << "\n";
	*/
	int x;
	cin >> x;

}