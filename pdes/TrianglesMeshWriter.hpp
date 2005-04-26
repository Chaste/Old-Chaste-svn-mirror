#ifndef _TRIANGLESMESHWRITER_HPP_
#define _TRIANGLESMESHWRITER_HPP_

#include "AbstractMeshWriter.hpp"

class TrianglesMeshWriter : public AbstractMeshWriter
{
public:
	TrianglesMeshWriter(std::string pathBaseName, int dimension);
	void WriteFiles();
	virtual ~TrianglesMeshWriter();
};

#endif //_TRIANGLESMESHWRITER_HPP_
