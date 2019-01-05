#include "MutationFinder.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
namespace MutationFinder{

	void output_to_file(std::string path, std::vector<MutationOutput> output){

		std::ofstream file(path);

		if (file.is_open()){
			for (auto line : output){

				char mut=' ';

				switch (line.mutation)
					{
					case INSERT:
						mut = 'I';
						break;
					case DELETE:
						mut = 'D';
						break;
					case SUBSTITUTION:
						mut = 'X';
						break;
					default:
						break;
					}

				file << mut << "," << line.position << "," << line.sign << "\n";

				}
			}

		file.close();

		}


	void map_mutations(int position, std::vector<std::string> aligned_strings, std::map<int, std::vector<MutationOutput>> &output){

		if (aligned_strings.size() == 0){
			return;
			}
		auto ref = aligned_strings[0];
		auto seq = aligned_strings[1];

		int insertions = 0;
		//checking for mutations
		for (int i = 0; i < ref.size(); i++){
				if (ref[i] == '-'){
					insertions++;
					auto pos = position + i - insertions;
					output[pos].push_back(MutationOutput(INSERT, seq[i], pos));
				}
				else if (seq[i] == '-'){
					auto pos = position + i - insertions;
					output[pos].push_back(MutationOutput(DELETE, '-', pos));
				}
			
				else if (ref[i] != seq[i]){
					auto pos = position + i - insertions;
					output[pos].push_back(MutationOutput(SUBSTITUTION, seq[i], pos));
				}
				else if (ref[i] == seq[i]){
					auto pos = position + i - insertions;
					output[pos].push_back(MutationOutput(NONE, seq[i], pos));
					}
		}
	}

	std::vector<MutationOutput> MapToVector(std::map<int, std::vector<MutationOutput>> outputMap){

		std::vector < MutationOutput > output;
		for (auto mutation : outputMap){

			std::unordered_map < MutationOutput, int > map;
			auto vm = mutation.second;
			for (auto mo : vm){
				map[mo]++;
				}
			/*for (auto mout : map){
				char mut = 'M';

				switch (mout.first.mutation)
					{
					case INSERT:
						mut = 'I';
						break;
					case DELETE:
						mut = 'D';
						break;
					case SUBSTITUTION:
						mut = 'X';
						break;
					default:
						break;
					}



				std::cerr << mut << " " << mout.first.position << " " << mout.first.sign << " count" << mout.second<<"\n";
				}*/
			MutationOutput most;
			int max=0;
			//finds most often mutations
			for (auto mut : map){
				if (mut.second > max){
					max = mut.second;
					most = mut.first;
				}
			}
			if (most.mutation != NONE){
				output.push_back(most);
				}

		}

		return output;
	}

	MutationOutput::MutationOutput(){
		sign = ' ';
		position = 0;
		mutation = INSERT;

		}

	MutationOutput::MutationOutput(Mutation mutatio, char sig, int positio){
		mutation = mutatio;
		sign = sig;
		position = positio;

		}
}