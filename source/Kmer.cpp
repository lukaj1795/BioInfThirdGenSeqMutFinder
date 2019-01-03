#include <string>
#include <map>
#include"Kmer.h"


Kmer::Kmer(std::string s, int p, int i) {
	identifier = i;
	string = s;
	position = p;
	ordering_number_for_string = 0;

	std::map<char, int> order_map_odd = {
	{ 'C', 0 },
	{ 'A', 1 },
	{ 'T', 2 },
	{ 'G', 3 }
	};
	std::map<char, int> order_map_even = {
	{ 'C', 3 },
	{ 'A', 2 },
	{ 'T', 1 },
	{ 'G', 0 }
	};

	/*We assign the values 0, 1, 2, 3 to C, A, T, G, respectively,
	for the odd numbered bases of k-mers, and reverse the
	ordering for even numbered bases*/
	/*if the base is even numbered*/
	int j = string.size()-1;
	for (std::string::size_type i = 0; i < string.size(); ++i) {
		if (i % 2 == 0) { //base is even numbered inside k-mer
			ordering_number_for_string += order_map_even.at(string[i])/ (std::pow(10,i));
		}
		else { //base is odd numbered inside k-mer
			ordering_number_for_string += order_map_odd.at(string[i])/ (std::pow(10,i));
		}
	}

}

Kmer::Kmer() {
	identifier = 0;
	position = 0;
	string = "";
}
bool Kmer::is_equal_to(Kmer kmer2) {
	if (this->string == kmer2.string && this->position == kmer2.position && this->identifier == kmer2.identifier) {
		return true;
	}
	return false;
}
