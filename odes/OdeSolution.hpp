#ifndef _ODESOLUTION_HPP_
#define _ODESOLUTION_HPP_


#include <vector>

class OdeSolution
{
	public:
	//OdeSolution();
	int mNumberOfTimeSteps;
	std::vector<double> mTime;
	std::vector<std::vector<double> > mSolutions;
};

#endif //_ODESOLUTION_HPP_
