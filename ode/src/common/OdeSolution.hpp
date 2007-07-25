/**
 * OdeSolution.  Sets us the class of ODE solutions, including a function that
 * allows us to save the output data to file.
 */
#ifndef _ODESOLUTION_HPP_
#define _ODESOLUTION_HPP_

#include <vector>
#include <fstream>
#include <cassert>
#include "ColumnDataWriter.hpp"
#include "AbstractOdeSystem.hpp"
#include <sstream>

class OdeSolution
{
private:
    unsigned mNumberOfTimeSteps;	/**< Variable for the number of timesteps */
    
    std::vector<double> mTimes; /**< A vector of times at each timestep. */
    std::vector<std::vector<double> > mSolutions;  /**< Solutions for each variable at each timestep. */
    
public:

    unsigned GetNumberOfTimeSteps(void)
    {
        return mNumberOfTimeSteps;
    }
    
    void SetNumberOfTimeSteps(unsigned num_timesteps)
    {
        mNumberOfTimeSteps = num_timesteps;
        mTimes.reserve(num_timesteps+1);
        mSolutions.reserve(num_timesteps+1);
    }
    
    std::vector<double> GetVariableAtIndex(unsigned index)
    {
        std::vector<double> answer;
        answer.reserve(mSolutions.size());
        for (unsigned i=0; i<mSolutions.size(); i++)
        {
            answer.push_back(mSolutions[i][index]);
        }
        return answer;
    }
    
    std::vector<double>& rGetTimes()
    {
        return mTimes;
    }
    
    std::vector<std::vector<double> >& rGetSolutions()
    {
        return mSolutions;
    }


    /**
     *  Write the data to a file. 
     *   @param pOdeSystem The ode system solved to obtain these results (needed for variable 
     *    names and units). 
     */ 
    void WriteToFile(std::string directoryName,
                     std::string baseResultsFilename,
                     AbstractOdeSystem* pOdeSystem,
                     std::string timeUnits,
                     unsigned stepPerRow = 1,
                     bool cleanDirectory = true)
    {
        assert(stepPerRow > 0);
        assert(mTimes.size()>0);
        assert(mTimes.size()==mSolutions.size());
     
        // Write data to a file using ColumnDataWriter
        ColumnDataWriter writer(directoryName,baseResultsFilename,cleanDirectory);

        int time_var_id = writer.DefineUnlimitedDimension("Time",timeUnits);
    
        std::vector<int> var_ids;
        for (unsigned i=0; i<pOdeSystem->rGetVariableNames().size(); i++)
        {
            var_ids.push_back(writer.DefineVariable(pOdeSystem->rGetVariableNames()[i],
                                                    pOdeSystem->rGetVariableUnits()[i]));
        }
        writer.EndDefineMode();

        
        for (unsigned i=0; i<mSolutions.size(); i+=stepPerRow)
        {
            writer.PutVariable(time_var_id, mTimes[i]);
            for (unsigned j=0; j<var_ids.size(); j++)
            {
                writer.PutVariable(var_ids[j], mSolutions[i][j]);
            }
            writer.AdvanceAlongUnlimitedDimension();
        }
        writer.Close();
    }
};

#endif //_ODESOLUTION_HPP_
