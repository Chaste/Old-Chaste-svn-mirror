/**
 * Concrete version of the AbstractMeshReader class.
 * A FemlabMeshReader takes the file names of a set of Femlab mesh files.
 * Once constructed the public methods of the AbstractMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */
#ifndef _FEMLABMESHREADER_H_
#define _FEMLABMESHREADER_H_

#include "AbstractMeshReader.cpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class FemlabMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>

{
private:
    std::vector<std::vector<double> > TokenizeStringsToDoubles(
        std::vector<std::string> rawData);
        
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(
        std::vector<std::string> rawData,
        unsigned dimensionOfObject);
public:
    FemlabMeshReader(std::string pathBaseName, std::string nodeFileName, std::string elementFileName, std::string edgeFileName);
    virtual ~FemlabMeshReader();
};

#endif //_FEMLABMESHREADER_H_
