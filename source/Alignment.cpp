#include "Alignment.h"
#include "Genome.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>

const int MATCH = 0;
const int INSERT = 1;
const int DELETE = 2;
const int window = 4;
namespace Alignment {

		int count = 0;
		std::vector<int> matrix;
		std::vector<int> operations;
		std::vector<int> operation{ 0, 0, 0 };

		int row;
		int column;
		int size;

		std::pair<int, std::vector<std::string>> Align(Kmer kmer_ref, Kmer kmer_seq) {


			std::vector<std::string> aligned_string;

			aligned_string.reserve(2);


			for (int i = 0; i < size; i++) {
				matrix[i] = 0;
				operations[i] = 0;
			}

			for (int i = 1; i < row; i++) {

				for (int j =std::max(1,i-window); j < column; j++) {


					if ((abs(i - j) > window)) {
						//std::cerr << (abs(i - j)) << "\n";
						//std::cerr << index << "\n";
						continue;
					}
					auto index = i * (column)+j;
					int match_award = -Genome::SUB+(kmer_seq.string[j - 1] == kmer_ref.string[i - 1])*(4+Genome::SUB);


					operation[MATCH] = matrix[index - column - 1] + match_award;
					operation[INSERT] = matrix[index - 1] - Genome::INSERT;
					operation[DELETE] = matrix[index - column] - Genome::DELETE;

					/*auto res=std::max(operation.begin(),operation.end());
					matrix[index] = *res;
					operations[index] = std::distance(operation.begin(),res);*/
					//std::cerr << matrix[index]<<"\n";
					//std::cerr << operations[index]<<"\n";
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

			if (max >= (row - 2) * 4 - Genome::DELETE) { //if only 1 mutation happens
				int index = bestIndex;
				std::string ref = "";
				std::string seq = "";

				///get first character
				ref += kmer_ref.string[index / column];
				seq += kmer_seq.string[index % column];
				count = 0;
				//backtracking
				Backtrack(operations, matrix, index, kmer_ref.string, kmer_seq.string, ref, seq);
				///we need to get the string in right order
				std::reverse(ref.begin(), ref.end());
				std::reverse(seq.begin(), seq.end());


				///we don't want to have excess characters
				int len = std::min(ref.length() - 1, seq.length() - 1);
				//std::cout << len << "\n";
				aligned_string.push_back(ref.substr(0, len));
				aligned_string.push_back(seq.substr(0, len));
				//std::cout <<"ref:	"<< kmer_ref.string << "\n";
				//std::cout <<"seq:	"<< kmer_seq.string << "\n" << "\nAlignment\n";
				//std::cout <<"ref:	"<< ref/*.substr(0, len) */ << '\n';
				//std::cout <<"seq:	"<< seq/*.substr(0, len)*/ << '\n' << "\n";
				//offset = index/column - count/column;
			}
			return std::move(std::make_pair(count / column, aligned_string));
		}


		void Backtrack(std::vector<int>operations, std::vector<int>matrix, int index, std::string const &ref, std::string const &seq, std::string &ref_align, std::string &seq_align) {

			while (matrix[index] != 0) {
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
			count = index;
		}


		void init_vectors(int k) {

			size = (k + 1)*(k + 1);
			row = k + 1;
			column = k + 1;
			matrix = std::vector<int>(size);
			operations = std::vector<int>(size);

		}


	}
