#include "OdeSolution.hpp"

void OdeSolution::SaveToFile(char * outputfile)
{
	std::ofstream file(outputfile);
	int numtimesteps = mTime.capacity();
	int numvariables= mSolutions[0].capacity();
	for(int i=0; i<numtimesteps; i++)
		{
			file << mTime[i] ;
			for(int k=0; k<numvariables; k++)
			{
				file  << "\t" << mSolutions[i][k];
			}
			file << std::endl;
		}
		file.close();
}


