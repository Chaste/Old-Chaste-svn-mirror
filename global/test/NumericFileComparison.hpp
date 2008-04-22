#ifndef NUMERICFILECOMPARISON_HPP_
#define NUMERICFILECOMPARISON_HPP_

#include "OutputFileHandler.hpp"

class NumericFileComparison
{ 
    /**
     * Compare files of numbers to see if they are to within a given tolerance. 
     */
private:
    std::ifstream file1, file2;
public:
    NumericFileComparison(std::string fileName1, std::string fileName2)
    {
        std::ifstream file1(fileName1.c_str(), std::ios::in);
        // If it doesn't exist - throw exception
        if (!file1.is_open())
        {
            EXCEPTION("Couldn't open info file: " + fileName1);
        }
        std::ifstream file2(fileName2.c_str(), std::ios::in);
        // If it doesn't exist - throw exception
        if (!file2.is_open())
        {
            EXCEPTION("Couldn't open info file: " + fileName2);
        }
    }
    bool CompareFiles( double absTolerance=DBL_EPSILON)
    {
        double data1, data2;
        unsigned failures=0;
        double max_error=0.0;
        unsigned max_failures=10;
        
        while (file1>>data1 && file2>>data2)
        {
            double error=fabs(data1 - data2);
            if ( error > absTolerance )
            {
                failures++;
                //Force CxxTest error
                TS_ASSERT_DELTA(data1, data2, absTolerance);
                if (error > max_error)
                {
                    max_error=error;
                }
            }
            if (failures > max_failures)
            {
                break;//Don't clog the screen
            }
        }
        //Can we read any more?
        if(file1>>data1 || file2>>data2)
        {
            EXCEPTION("Files have different lengths");
        }
        //Force CxxTest error if there were any major differences
        TS_ASSERT_LESS_THAN(max_error, absTolerance);
        
        return (failures==0);   
    }                       
};
#endif /*NUMERICFILECOMPARISON_HPP_*/
