#ifndef _MESHALYZERMESHWRITER_HPP_
#define _MESHALYZERMESHWRITER_HPP_

#include "AbstractMeshWriter.hpp"


class MeshalyzerMeshWriter : public AbstractMeshWriter
{
public:
	MeshalyzerMeshWriter(std::string pathBaseName, bool setCoolGraphics=false);
	void WriteFiles();
	virtual ~MeshalyzerMeshWriter();

};

#endif //_MESHALYZERMESHWRITER_HPP_
