#include"FASTAUtility.h"
#include<iostream>
namespace FASTAUtility{
	std::vector<std::string> ReadFasta(std::ifstream &file){
		std::string line, sequence;
		auto i = 0;
		std::vector<std::string> vector;
		while (std::getline(file, line)){
			//sanity check
			if (line.empty()){
				continue;
			}
			//new sequence
			if (line[0] == '>'){
				
				//we don't want an empty sequence as the first one
				if (i != 0){
					vector.push_back(sequence);
				}
				sequence = "";
				i++;
			}
			//append to sequence
			else{
				sequence += line;

			}
		}
		vector.push_back(sequence);
		return vector;
	}
}