#include<string>

class Genome{

public:

	static const int SUB=2;
	static const int INSERT=2;
	static const int DELETE=3;
	int identifier; /*identifier*/

	std::string genomeString;

	Genome(int i);

	
};