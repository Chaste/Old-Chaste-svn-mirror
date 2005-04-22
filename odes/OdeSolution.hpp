/**
 * OdeSolution.  Sets us the class of ODE solutions, including a function that 
 * allows us to save the output data to file.
 */
#ifndef _ODESOLUTION_HPP_
#define _ODESOLUTION_HPP_

#include <vector>
#include <fstream>

class OdeSolution
{
	public:
	int mNumberOfTimeSteps;	/** Variable for the number of timesteps */
	std::vector<double> mTime; /** Sets up a vector of time. */
	std::vector<std::vector<double> > mSolutions;  /** Sets up the solutions. */
	void SaveToFile(char * outputfile);
};

#endif //_ODESOLUTION_HPP_
