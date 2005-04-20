#ifndef _ABSTRACTMESHREADER_HPP_
#define _ABSTRACTMESHREADER_HPP_

#include <string>

class AbstractMeshReader
{
	private:
		int mNumElements;
		bool mSuccess;
		std::string mPathBaseName;
	public:
	bool IsReaderSuccess(){return(mSuccess);}
	AbstractMeshReader(std::string pathBaseName)
	{
		mPathBaseName=pathBaseName;
		//This is where we do something
		mSuccess=true;
	}
};
#endif //_ABSTRACTMESHREADER_HPP_
