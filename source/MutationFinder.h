#include<vector>
#include<map>
//namespace used for finding mutations
namespace MutationFinder{

	enum Mutation {
		INSERT,
		DELETE,
		SUBSTITUTION,
		NONE
		};
	//class that represents mutation output
	class MutationOutput
		{
		public:
			MutationOutput();
			Mutation mutation; //mutation that happened
			char sign; // sign that was inserted, matched, deleted or substitution sign
			int position; // position in which it happened
			MutationOutput(Mutation, char, int);

			bool operator==(const MutationOutput &other) const
				{
				return (mutation == other.mutation && sign==other.sign && position==other.position );
				}

			friend struct std::hash<MutationOutput>;
		
		};


	//maps mutations to output map
	void map_mutations(int position, std::vector<std::string> aligned_strings,std::map<int,std::vector<MutationOutput>> &output);
	//outputs given vector to file in given path
	void output_to_file(std::string path, std::vector<MutationOutput> output);
	//creates vector with most probable mutations
	std::vector<MutationOutput> MapToVector(std::map<int, std::vector<MutationOutput>> outputMap);
	
}


///add hash template to std for mutation output, needed for unoredered_map of outputs
namespace std
	{
	template <>
	struct hash<MutationFinder::MutationOutput>
		{
		size_t operator()(const MutationFinder::MutationOutput& k) const{
			return ((hash<int>()(k.sign) ^ (hash<char>()(k.sign) << 1)^hash<int>()(k.mutation)) >> 1);
			}
		};			

	}