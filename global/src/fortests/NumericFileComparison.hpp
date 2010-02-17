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
#define A_WORD DBL_MAX
#define NOTHING_TO_READ DBL_MIN
/**
 * Compare files of numbers to see if they match to within a given tolerance.
 */
class NumericFileComparison
{
private:
    std::string mFilename1; /**< First filename */
    std::string mFilename2; /**< Second filename */

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
    NumericFileComparison(std::string fileName1, std::string fileName2):
        mFilename1(fileName1),
        mFilename2(fileName2)
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
     * @param ignoreFirstFewLines  How many lines to ignore from the comparison
     * 
     */
    bool CompareFiles(double absTolerance=DBL_EPSILON, unsigned ignoreFirstFewLines=0)
    {
        double data1;
        double data2;
        unsigned failures = 0;
        double max_error = 0.0;
        unsigned max_failures = 10;
        
        for (unsigned line_number=0; line_number<ignoreFirstFewLines; line_number++)
        {
            char buffer[256];
            mpFile1->getline(buffer, 256);
            mpFile2->getline(buffer, 256);
            TS_ASSERT(!mpFile1->fail()); //Here we are assuming that there a least "ignoreFirstFewLines" lines
            TS_ASSERT(!mpFile2->fail()); // and that they are lines of no more than 256 characters
        }
        
        do 
        {
            if (!(*mpFile1>>data1))
            {
                //Cannot read the next token from file as a number, so try a word instead
                std::string word;
                mpFile1->clear();//reset the "failbit"
                if (*mpFile1>>word)
                {
                    data1=A_WORD;
                }
                else
                {
                    mpFile1->clear();//reset the "failbit"
                    data1=NOTHING_TO_READ;          
                }
            }
            if (!(*mpFile2>>data2))
            {
                //Cannot read the next token from file as a number, so try a word instead
                std::string word;
                mpFile2->clear();//reset the "failbit"
                if (*mpFile2>>word)
                {
                    data2=A_WORD;
                }
                else
                {
                    mpFile2->clear();//reset the "failbit"
                    data2=NOTHING_TO_READ;          
                }
            }
            
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
        while (data1 != NOTHING_TO_READ && data2 != NOTHING_TO_READ); //If either is a NOTHING_TO_READ, then it means that there's nothing to read from the file

        // Force CxxTest error if there were any major differences
        TS_ASSERT_LESS_THAN(max_error, absTolerance);
        //If that assertion tripped...
        if (max_error >= absTolerance)
        {
#define COVERAGE_IGNORE            
            //Report the paths to the files
            TS_TRACE("Files " + mFilename1 + " and " + mFilename2 + " numerically differ.");
#undef COVERAGE_IGNORE            
        }
        return (failures==0);
    }
};

#endif /*NUMERICFILECOMPARISON_HPP_*/
