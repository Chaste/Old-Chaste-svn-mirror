#ifndef _ODESOLUTION_HPP_
#define _ODESOLUTION_HPP_

#include <vector>
#include <fstream>

class OdeSolution
{
	public:
	int mNumberOfTimeSteps;
	std::vector<double> mTime;
	std::vector<std::vector<double> > mSolutions;
	void SaveToFile(char * outputfile);
};

#endif //_ODESOLUTION_HPP_
