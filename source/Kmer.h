#ifndef KMER_H
#define KMER_H
#endif // !KMER_H

#include <map>
#include <string>

class Kmer { /*(i, p) identifying the sequence i and the position p within sequence at
which s appears. (s, i, p) is called a k-mer triple */

public:
	std::string string; /*actual k-mer*/
	
	int position; /*position inside the string*/
	int identifier; /*sequence i*/

	int ordering_number_for_string;/*We assign the values 0, 1, 2, 3 to C, A, T, G, respectively,
	for the odd numbered bases of k-mers, and reverse the
	ordering for even numbered bases*/
	bool is_equal_to(Kmer kmer2);
	Kmer(std::string string, int p, int i);
	//Kmer();
	//Kmer(Kmer&& other);

	bool operator<(const Kmer& rhs) const
	{
		return ordering_number_for_string < rhs.ordering_number_for_string;
	}


	/*Kmer& operator=(Kmer&& other){

		if (this != &other){
			

			
			
			}

		return *this;
		}*/

};
