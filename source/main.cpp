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
			cout <<"number of references:	"<< references.size() << "\n";

		}
		if (file2.is_open()){
			V = ReadFasta(file2);
			cout <<"number of sequences:	"<< V.size() << "\n";

			}
		file.close();
		file2.close();
	}
	auto finish = chrono::high_resolution_clock::now();
	cout <<"Reading input files:	"<< std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n\n";
	//cout << reference << "\n";
	//cout << V[0] << "\n";


	Genome g1= Genome(1);
	Genome g2 = Genome(2);
	//g1.genomeString = "GTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGATGCATCAGCGAGCAGCGATCGTAAGACCGCTAGCTAAGCTCGGCCTACGATAACGACATCAGCTACGATGCATCGATCTGATCGAGCATGCTGAGCAGCGTACTATGCGTAGTCATGCTGAGTGTCTTGGTCAGCAAAATGCATCGATCGACATGGTGTTCGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGAGGTCGTCACGTTACTCACAAGCAT";
	//g2.genomeString = "CGATCGTAAGACCGCTAGCTAAGCTCGGGGCATCGACTACAAAATTTCGAGTGATCGATGCCAGA";
	start = chrono::high_resolution_clock::now();
	g1.genomeString = references[0];
	g2.genomeString = V[0];
	
	Kmer_extraction kmer = Kmer_extraction(20,20);
	auto kmer1=kmer.extract(&g1);
	finish = chrono::high_resolution_clock::now();
	cout << "k_mer extraction reference:	"<< std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";

	start = chrono::high_resolution_clock::now();
	std::map<int,int> pokrivenost1;
	int counter = 1;
	for (auto i : V) {
		cout << "\n\nstring:	" << counter << "	duljina:	" << i.length();
		g2.genomeString = i;
		auto kmer2 = kmer.extract(&g2);
		int position = mapping::map_to_reference(&kmer1, &kmer2);
		if (position != -99999) {
			cout << "\npozicija unutar reference gdje bi se trebali poceti preklapati:	" << position << "	+-10 mjesta" << "\n";
			pokrivenost1.insert(std::pair<int,int>(counter,position));
		}
		counter++;
	}
	finish = chrono::high_resolution_clock::now();
	cout <<"Alignment:	"<< std::chrono::duration_cast<chrono::nanoseconds>(finish - start).count() / 1e6 << " ms\n";
	
	int counter = 1;
	for (auto i : pokrivenost1) {
		cout << "\nstring:	" <<i.first<<"	position:	"<< i.second;
	}
	/*
	cout << "\n\n" << "string 1: " << g1.genomeString << "\n\n";
	for (auto i : kmer1)
		std::cout <<"kmer: "<< i.string << " vrijednost: " << i.ordering_number_for_string<<" pozicija: " << i.position << " sekvenca: " << i.identifier <<"\n";
	
	cout <<"\n\n"<< "string 2: " << g2.genomeString << "\n\n";
	for (auto i : kmer2)
		std::cout << "kmer: " << i.string << " vrijednost: " << i.ordering_number_for_string << " pozicija: " << i.position << " sekvenca: " << i.identifier << "\n";
	*/
	
	
	int x;
	cin >> x;

}