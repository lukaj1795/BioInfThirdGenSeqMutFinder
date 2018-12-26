#include"Genome.h"
#include"Kmer_extraction.h"
#include"Kmer.h"
#include "FASTAUtility.h"
#include <iostream>
#include <string>
#include <vector>
#include<fstream>

using namespace std;
using namespace FASTAUtility;

void main(int argc,char** argv){

	//reading fasta files
	vector<string> V;
	if (argc == 2){
		ifstream file(argv[1]);
		if (file.is_open()){
			V = ReadFasta(file);
			cout << V.size() << "\n";

		}
		
	}
	
	//cout << V[0] << "\n";
	Genome g1= Genome(1);
	Genome g2 = Genome(2);
	Genome g3 = Genome(3);
	g3.genomeString = "AGGCTGTCACGTCAGTGTGCCAACGTGCAACGCCTGACCTGACTGGGTCACTGACTG";
	g2.genomeString = "CGTCAGTGTGCCAACGTGCAACGCCT";
	g1.genomeString = "AACGCCTGACCTGACTGGGTCA";
	//g1.genomeString = V[0];
	Kmer_extraction kmer = Kmer_extraction(4,5);
	std::vector<Kmer> svi_kmerovi_sekvence;

	cout << "string 1: " << g1.genomeString << "\n\n";
	svi_kmerovi_sekvence = kmer.extract(&g1);
	for (auto i : svi_kmerovi_sekvence)
		std::cout <<"kmer: "<< i.string << " vrijednost: " << i.ordering_number_for_string<<" pozicija: " << i.position << " sekvenca: " << i.identifier <<"\n";

	cout <<"\n\n"<< "string 2: " << g2.genomeString << "\n\n";
	svi_kmerovi_sekvence = kmer.extract(&g2);
	for (auto i : svi_kmerovi_sekvence)
		std::cout << "kmer: " << i.string << " vrijednost: " << i.ordering_number_for_string << " pozicija: " << i.position << " sekvenca: " << i.identifier << "\n";

	cout <<"\n\n"<< "string 3: " << g3.genomeString << "\n\n";
	svi_kmerovi_sekvence = kmer.extract(&g3);
	for (auto i : svi_kmerovi_sekvence)
		std::cout << "kmer: " << i.string << " vrijednost: " << i.ordering_number_for_string << " pozicija: " << i.position << " sekvenca: " << i.identifier << "\n";

	int x;
	cin >> x;

}