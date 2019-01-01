#include "Kmer_mapping.h"

Kmer_mapping::Kmer_mapping(int x, int y) {
	sequence_position = x;
	reference_positions.push_back(y);
}

void Kmer_mapping::add_position(int y) {
	reference_positions.push_back(y);
}

