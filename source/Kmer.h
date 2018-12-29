#ifndef KMER_H
#define KMER_H
#endif // !KMER_H


#include <string>

class Kmer { /*(i, p) identifying the sequence i and the position p within sequence at
which s appears. (s, i, p) is called a k-mer triple */

public:
	std::string string; /*actual k-mer*/
	
	int position; /*position inside the string*/
	int identifier; /*sequence i*/
	
	std::string ordering_number_for_string;/*We assign the values 0, 1, 2, 3 to C, A, T, G, respectively,
	for the odd numbered bases of k-mers, and reverse the
	ordering for even numbered bases*/

	Kmer(std::string string, int p, int i);
};
