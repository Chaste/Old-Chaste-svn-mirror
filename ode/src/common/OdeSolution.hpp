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
	int mNumberOfTimeSteps;	/**< Variable for the number of timesteps */
	
    std::vector<double> mTimes; /**< A vector of times at each timestep. */
    std::vector<std::vector<double> > mSolutions;  /**< Solutions for each variable at each timestep. */
 
public:

	int GetNumberOfTimeSteps(void)
	{
		return mNumberOfTimeSteps;
	}
	void SetNumberOfTimeSteps(int num_timesteps)
	{
		mNumberOfTimeSteps = num_timesteps;
		mTimes.reserve(num_timesteps+1);
		mSolutions.reserve(num_timesteps+1);
	} 
    
    std::vector<double> GetVariableAtIndex(int index)
    {
        std::vector<double> answer;
        for (unsigned i=0; i<mSolutions.size(); i++){
            answer.push_back(mSolutions[i][index]);
        }
        return(answer);
    }
    std::vector<double>& rGetTimes()
    {
        return mTimes;
    }
	
    std::vector<std::vector<double> >& rGetSolutions()
    {
        return mSolutions;
    }
};

#endif //_ODESOLUTION_HPP_
