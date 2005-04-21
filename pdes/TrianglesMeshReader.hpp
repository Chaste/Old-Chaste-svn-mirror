#ifndef _TRIANGLESMESHREADER_H_
#define _TRIANGLESMESHREADER_H_

#include "AbstractMeshReader.hpp"
#include "Exception.hpp"

class TrianglesMeshReader : public AbstractMeshReader
{
private:
	std::vector<std::vector<double> > TokenizeStringsToDoubles(
											std::vector<std::string> rawData);
	std::vector<std::vector<int> > TokenizeStringsToInts(
											std::vector<std::string> rawData,
											int dimensionOfObject);											
public:
	TrianglesMeshReader(std::string pathBaseName);
	virtual ~TrianglesMeshReader();
};

#endif //_TRIANGLESMESHREADER_H_
