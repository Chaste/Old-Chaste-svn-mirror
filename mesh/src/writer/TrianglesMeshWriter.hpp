#ifndef _TRIANGLESMESHWRITER_HPP_
#define _TRIANGLESMESHWRITER_HPP_

#include "AbstractMeshWriter.hpp"
#include "OutputFileHandler.hpp"

class TrianglesMeshWriter : public AbstractMeshWriter
{
public:
	TrianglesMeshWriter(const std::string &rDirectory, const std::string &rBbaseName, const int &rDimension);
	void WriteFiles();
	virtual ~TrianglesMeshWriter();
};

#endif //_TRIANGLESMESHWRITER_HPP_
