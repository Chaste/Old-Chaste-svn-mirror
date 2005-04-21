#include "AbstractMeshReader.hpp"


std::vector<std::string> AbstractMeshReader::GetRawDataFromFile(std::string fileName)
{
	std::vector<std::string> RawDataFromFile;

	std::ifstream dataFile(fileName.c_str());
	
	if (!dataFile.is_open())
	{
		throw Exception("Could not open data file "+fileName+" .");
	}
	
	std::string RawLineFromFile;
	getline(dataFile, RawLineFromFile);
	
	while(dataFile){
		int hashLocation=RawLineFromFile.find('#',0);
		if (hashLocation >= 0)
		{
			RawLineFromFile=RawLineFromFile.substr(0,hashLocation);
		}
		int notBlankLocation=RawLineFromFile.find_first_not_of(" \t",0);
		if (notBlankLocation >= 0)
		{
			RawDataFromFile.push_back(RawLineFromFile);	
		}
		getline(dataFile, RawLineFromFile);
	}


	dataFile.close();
	
	return(RawDataFromFile);
}
