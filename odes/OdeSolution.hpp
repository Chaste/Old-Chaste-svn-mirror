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
	private:
	int mNumberOfTimeSteps;	/** Variable for the number of timesteps */
	
	public:
	int GetNumberOfTimeSteps(void) {return mNumberOfTimeSteps;}
	void SetNumberOfTimeSteps(int num_timesteps) {mNumberOfTimeSteps = num_timesteps;} 
	std::vector<double> mTime; /** Sets up a vector of time. */
	std::vector<std::vector<double> > mSolutions;  /** Sets up the solutions. */
	
};

#endif //_ODESOLUTION_HPP_
