#include<vector>
#include<map>

namespace MutationFinder{

	enum Mutation {
		INSERT,
		DELETE,
		SUBSTITUTION
		};

	class MutationOutput
		{
		public:
			MutationOutput();
			Mutation mutation;
			char sign;
			int position;
			MutationOutput(Mutation, char, int);

			bool operator==(const MutationOutput &other) const
				{
				return (mutation == other.mutation && sign==other.sign && position==other.position );
				}

			friend struct std::hash<MutationOutput>;
		
		};


	//maps mutations to output map
	void map_mutations(int position, std::vector<std::string> aligned_strings,std::map<int,std::vector<MutationOutput>> &output);
	//output given vector to file in given path
	void output_to_file(std::string path, std::vector<MutationOutput> output);
	//creates vector with most probalble mutations
	std::vector<MutationOutput> MapToVector(std::map<int, std::vector<MutationOutput>> outputMap);
	
}


///add hash template to std for mutation output, needed for unoredered_map of outputs
namespace std
	{
	template <>
	struct hash<MutationFinder::MutationOutput>
		{
		size_t operator()(const MutationFinder::MutationOutput& k) const{
			// Compute individual hash values for two data members and combine them using XOR and bit shifting
			return ((hash<int>()(k.sign) ^ (hash<char>()(k.sign) << 1)^hash<int>()(k.mutation)) >> 1);
			}
		};			

	}