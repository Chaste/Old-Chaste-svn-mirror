#ifndef _TRIANGLESMESHWRITER_HPP_
#define _TRIANGLESMESHWRITER_HPP_

#include "AbstractMeshWriter.cpp"
#include "OutputFileHandler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
public:
    TrianglesMeshWriter(const std::string &rDirectory,
                        const std::string &rBbaseName,
                        const bool clearOutputDir=true);
    void WriteFiles();
    void WriteElementsAsFaces();
    void WriteFacesAsEdges();
    virtual ~TrianglesMeshWriter();
};

#endif //_TRIANGLESMESHWRITER_HPP_
