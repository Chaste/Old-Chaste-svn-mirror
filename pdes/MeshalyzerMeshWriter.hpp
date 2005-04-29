#ifndef _MESHALYZERMESHWRITER_HPP_
#define _MESHALYZERMESHWRITER_HPP_

#include "AbstractMeshWriter.hpp"


class MeshalyzerMeshWriter : public AbstractMeshWriter
{
public:
	MeshalyzerMeshWriter(std::string pathBaseName);
	void WriteFiles();
	void SetCoolGraphicsFormat();
	virtual ~MeshalyzerMeshWriter();

};

#endif //_MESHALYZERMESHWRITER_HPP_
