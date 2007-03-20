#ifndef _MESHALYZERMESHWRITER_HPP_
#define _MESHALYZERMESHWRITER_HPP_

#include "AbstractMeshWriter.cpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MeshalyzerMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
public:
    MeshalyzerMeshWriter(const std::string &rDirectory,
                         const std::string &rBaseName,
                         const bool &rCleanDirectory=true,
                         const bool &rSetCoolGraphics=false);
    void WriteFiles();
    virtual ~MeshalyzerMeshWriter();
    
};

#endif //_MESHALYZERMESHWRITER_HPP_
