/** 
 * Concrete version of the AbstractMeshReader class.
 * A TrianglesMeshReader takes the base name of a set of Triangles or
 * Tetgen mesh files (ie. the path and name of the files without the suffices).
 * Once constructed the public methods of the AbstractMeshReader 
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */
#ifndef _TRIANGLESMESHREADER_H_
#define _TRIANGLESMESHREADER_H_

#include "AbstractMeshReader.hpp"
#include "global/src/Exception.hpp"

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
