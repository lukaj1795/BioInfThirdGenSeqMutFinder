#include<string>

//class that represent a genome
class Genome{

public:

	static const int SUB=2;//cost for substitution
	static const int INSERT=3;//cost for insert
	static const int DELETE=3;//cost for delete
	int identifier; /*identifier*/

	std::string genomeString; //string in genome

	Genome(int i);

	
};