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
	if (argc == 2){
		ifstream file(argv[1]);
		if (file.is_open()){
			//V = ReadFasta(file);
			//cout << V.size() << "\n";

		}
		
	}
	auto finish = chrono::high_resolution_clock::now();
	//cout << std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";
	//cout << V[0] << "\n";
	Genome g1= Genome(1);
	Genome g2 = Genome(2);
	Genome g3 = Genome(3);
	Genome g4 = Genome(4);
	g3.genomeString = "AGGCTGTCACGTCAGTGTGCCAACGTGCAACGCCTGACCTGACTGGGTCACTGACTG";
	g2.genomeString = "CGTCAGTGTGCCAACGTGCAACGCCT";
	g1.genomeString = "AACGCCTGACCTGACTGGGTCA";
	g4.genomeString = "AACGCCTAGCCTGACTGGGTCA";
	//g1.genomeString = V[0];
	Kmer_extraction kmer = Kmer_extraction(4,5);
	Kmer_extraction k2 = Kmer_extraction(4,6);
	//std::vector<Kmer> svi_kmerovi_sekvence;
	auto kmer1=k2.extract(&g1);
	auto kmer4=kmer.extract(&g4);

	Alignment::Align(kmer1,kmer4);
	//cout << "string 1: " << g1.genomeString << "\n\n";
	start = chrono::high_resolution_clock::now();
	auto svi_kmerovi_sekvence = kmer.extract(&g1);
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