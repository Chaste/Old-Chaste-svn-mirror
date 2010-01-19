/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef NUMERICFILECOMPARISON_HPP_
#define NUMERICFILECOMPARISON_HPP_

#include "OutputFileHandler.hpp"

/**
 * Compare files of numbers to see if they match to within a given tolerance.
 */
class NumericFileComparison
{
private:
    std::ifstream* mpFile1; /**< First file */
    std::ifstream* mpFile2; /**< Second file */
public:
    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     * 
     * @param fileName1  first file
     * @param fileName2  second file
     */
    NumericFileComparison(std::string fileName1, std::string fileName2)
    {
        mpFile1 = new std::ifstream(fileName1.c_str());
        // If it doesn't exist - throw exception
        if (!mpFile1->is_open())
        {
            mpFile1 = NULL;
            EXCEPTION("Couldn't open file: " + fileName1);
        }

        mpFile2 = new std::ifstream(fileName2.c_str());
        // If it doesn't exist - throw exception
        if (!mpFile2->is_open())
        {
            mpFile1->close();
            mpFile1 = NULL;
            mpFile2 = NULL;
            EXCEPTION("Couldn't open file: " + fileName2);
        }
    }
    /**
     * Close the files being compared.
     */
    ~NumericFileComparison()
    {
        if (mpFile1)
        {
            mpFile1->close();
        }
        if (mpFile2)
        {
            mpFile2->close();
        }
        delete mpFile1;
        delete mpFile2;
    }
    /**
     * Compare the files.
     * 
     * @param absTolerance  absolute tolerance on difference between numbers.
     */
    bool CompareFiles(double absTolerance=DBL_EPSILON)
    {
        double data1, data2;
        unsigned failures = 0;
        double max_error = 0.0;
        unsigned max_failures = 10;
        bool empty_files = true;
        
        while (*mpFile1>>data1 && *mpFile2>>data2)
        {
            empty_files = false;
            double error = fabs(data1 - data2);
            if ( error > absTolerance )
            {
                failures++;
                // Force CxxTest error
                TS_ASSERT_DELTA(data1, data2, absTolerance);
                if (error > max_error)
                {
                    max_error = error;
                }
            }
            if (failures > max_failures)
            {
                break; // Don't clog the screen
            }
        }
        // Can we read any more?
        if (*mpFile1>>data1 || *mpFile2>>data2)
        {
            EXCEPTION("Files have different lengths");
        }
        // Force CxxTest error if there were any major differences
        TS_ASSERT_LESS_THAN(max_error, absTolerance);
        TS_ASSERT(!empty_files);
        return (failures==0 && !empty_files);
    }
};

#endif /*NUMERICFILECOMPARISON_HPP_*/
