#include "TrianglesMeshWriter.hpp"

TrianglesMeshWriter::TrianglesMeshWriter(std::string pathBaseName, int dimension)
{
	
	//Copy path and base name of files to private data
	mPathBaseName=pathBaseName;
	mDimension=dimension;
}

void
TrianglesMeshWriter::WriteFiles()
{
	
	
}


TrianglesMeshWriter::~TrianglesMeshWriter()
{
}
