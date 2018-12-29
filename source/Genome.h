#include<string>

class Genome{

public:

	static const int SUB=1;
	static const int INSERT=2;
	static const int DELETE=2;
	int identifier; /*identifier*/

	std::string genomeString;

	Genome(int i);

	
};