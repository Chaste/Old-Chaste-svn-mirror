#ifndef _TRIANGLESMESHREADER_H_
#define _TRIANGLESMESHREADER_H_

#include "AbstractMeshReader.hpp"
#include "Exception.hpp"

class TrianglesMeshReader : public AbstractMeshReader
{
public:
	TrianglesMeshReader(std::string pathBaseName);
	virtual ~TrianglesMeshReader();
};

#endif //_TRIANGLESMESHREADER_H_
