#include"Genome.h"
#include<iostream>
using namespace std;
void main(){
	Genome g=Genome();
	g.genomeString = "AACGCCTGA";
	cout << g.genomeString;
	cout << g.INSERT;
	int x;
	cin >> x;
}