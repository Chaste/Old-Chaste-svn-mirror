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

#include "AbstractMeshReader.cpp"
#include "Exception.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:
    std::vector<std::vector<double> > TokenizeStringsToDoubles(
        std::vector<std::string> rawData);
        
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(
        std::vector<std::string> rawData,
        unsigned dimensionOfObject);
    void ReadFacesAsElements(std::string pathBaseName);
    void ReadEdgesAsFaces(std::string pathBaseName);
    
public:
    TrianglesMeshReader(std::string pathBaseName);
    
    virtual ~TrianglesMeshReader();
};

#endif //_TRIANGLESMESHREADER_H_
