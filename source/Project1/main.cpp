#include"Genome.h"
#include"Kmer_extraction.h"
#include<iostream>
using namespace std;
void main(){
	Genome g=Genome();
	g.genomeString = "AACGCCTGACCTGACTGGGTCA";
	cout << g.genomeString << ' ';
	cout << g.INSERT << "\n";

	Kmer_extraction kmer = Kmer_extraction(2,4);
	
	std::vector<std::string> novo_polje = kmer.extract(g.genomeString);

	for (auto i : novo_polje)
		std::cout << i << ' ';

	int x;
	cin >> x;

}