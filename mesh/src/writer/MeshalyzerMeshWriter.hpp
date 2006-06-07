#ifndef _MESHALYZERMESHWRITER_HPP_
#define _MESHALYZERMESHWRITER_HPP_

#include "AbstractMeshWriter.hpp"


class MeshalyzerMeshWriter : public AbstractMeshWriter
{
public:
	MeshalyzerMeshWriter(const std::string &rDirectory, 
                         const std::string &rBaseName, 
                         const bool &rSetCoolGraphics=false);
	void WriteFiles();
	virtual ~MeshalyzerMeshWriter();

};

#endif //_MESHALYZERMESHWRITER_HPP_
