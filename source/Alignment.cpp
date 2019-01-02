#include "Alignment.h"
#include "Genome.h"
#include <iostream>
#include <algorithm>
#include <vector>

const int MATCH = 0;
const int INSERT = 1;
const int DELETE = 2;

namespace Alignment{
	std::vector<std::string> Align(Kmer reference, Kmer sequence){

		
		std::vector<std::string> aligned_string;
		auto kmer_ref = reference;
		auto kmer_seq = sequence;
		int n = kmer_ref.string.length();
		int m = kmer_seq.string.length();
				const int size = (n + 1)*(m + 1);
				auto row = n + 1;
				auto matrix = std::vector<int>(size);
				auto operations = std::vector<int>(size);
				int operation[3] = { 0 };
				auto column = m + 1;

				for (int i = 0; i < size; i++) {
					matrix[i] = 0;
					operations[i] = 0;
				}

				for (int i = 1; i < row; i++) {
					for (int j = 1; j < column; j++) {

						auto index = i * (column)+j;

						int match_award = -Genome::SUB;

						if (kmer_seq.string[j - 1] == kmer_ref.string[i - 1]) {
							match_award = 4;
						}
						operation[MATCH] = matrix[index - column - 1] + match_award;
						operation[INSERT] = matrix[index - 1] - Genome::INSERT;
						operation[DELETE] = matrix[index - column] - Genome::DELETE;

						for (int k = MATCH; k <= DELETE; k++) {

							if (operation[k] > matrix[index]) {
								matrix[index] = operation[k];
								operations[index] = k;
							}
						}
						//std::cout<<i<<" "<<j<<" "<< matrix[index] << "\t";

					}

				}

				///finding the highest score for local alignment
				int max = -1;
				int bestIndex = -1;
				for (int index = 0; index < size; index++) {
					if (matrix[index] > max) {
						bestIndex = index;
						max = matrix[index];
					}
				}

				if (max >= (n * 4 - 4 - 2 * 2+1)) { //if only 1 mutation happens
					int index = bestIndex;
					std::string ref = "";
					std::string seq = "";

					///get first character
					ref += kmer_ref.string[index / column];
					seq += kmer_seq.string[index % column];

					//backtracking
					Backtrack(operations, matrix, column, index, kmer_ref.string, kmer_seq.string, ref, seq);
					///we need to get the string in right order
					std::reverse(ref.begin(), ref.end());
					std::reverse(seq.begin(), seq.end());

					
					///we don't want to have excess characters
					int len = std::min(ref.length() - 1, seq.length() - 1);
					//std::cout << len << "\n";
					aligned_string.push_back(ref.substr(0,len));
					aligned_string.push_back(seq.substr(0,len));
					//std::cout <<"ref:	"<< kmer_ref.string << "\n";
					//std::cout <<"seq:	"<< kmer_seq.string << "\n" << "\nAlignment\n";
					//std::cout <<"ref:	"<< ref/*.substr(0, len) */ << '\n';
					//std::cout <<"seq:	"<< seq/*.substr(0, len)*/ << '\n' << "\n";
					
				}
		return aligned_string;
	} 


	void Backtrack(std::vector<int>operations, std::vector<int>matrix, int column, int index, std::string const &ref, std::string const &seq, std::string &ref_align, std::string &seq_align){
		while (matrix[index] != 0){
			switch (operations[index])
				{
				case MATCH:
					index = index - column - 1;
					ref_align += ref[index / column];
					seq_align += seq[index%column];
					break;
				case INSERT:
					index = index - 1;
					ref_align += '-';
					seq_align += seq[index%column];
					break;
				case DELETE:
					index = index - column;
					ref_align += ref[index / column];
					seq_align += '-';
					break;
				default:
					break;

				}
			}
		
		}
}