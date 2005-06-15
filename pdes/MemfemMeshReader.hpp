#ifndef _MEMFEMMESHREADER_HPP_
#define _MEMFEMMESHREADER_HPP_

/** 
 * Concrete version of the AbstractMeshReader class.
 * A MemfemMeshReader takes the base name of a set of Memfem 
 * mesh files (ie. the path and name of the files without the suffices).
 * Once constructed the public methods of the AbstractMeshReader 
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */

#include "AbstractMeshReader.hpp"
#include "common/Exception.hpp"

class MemfemMeshReader : public AbstractMeshReader
{
	
	private:
	std::vector<std::vector<double> > TokenizeStringsToDoubles(
											std::vector<std::string> rawData);
											
	std::vector<std::vector<int> > TokenizeStringsToInts(
											std::vector<std::string> rawData,
											int dimensionOfObject,
											bool readHeader);											

	public:
		MemfemMeshReader(std::string pathBaseName);
		virtual ~MemfemMeshReader();
};

#endif //_MEMFEMMESHREADER_HPP_
