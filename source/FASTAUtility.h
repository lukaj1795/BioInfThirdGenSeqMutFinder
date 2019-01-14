#include<string>
#include<vector>
#include<fstream>
//namespace used for providing some FASTA utility
namespace FASTAUtility{
//Used to read fasta file, returns sequences in order of appearance
 std::vector<std::string> ReadFasta(std::ifstream & file); 



};