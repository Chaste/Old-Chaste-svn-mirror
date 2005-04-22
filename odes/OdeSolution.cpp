/**
 * OdeSolution.  Sets us the function that prints the output data to a file.
 */
#include "OdeSolution.hpp"

/**
 * Solves a system of ODEs using the Forward Euler method
 * 
 * @param outputfile points to the output file* 
 *  
 * To be used in the form:
 * 
 *  solution.SaveToFile("Name_of_file");
*/

void OdeSolution::SaveToFile(char * pOutputfile)
{
	std::ofstream file(pOutputfile);
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


