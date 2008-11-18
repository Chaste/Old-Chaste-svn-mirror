/*

Copyright (C) University of Oxford, 2008

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
#ifndef PROGRESSREPORTER_HPP_
#define PROGRESSREPORTER_HPP_

#include <cassert>
#include "OutputFileHandler.hpp"


/**
 *  This class creates a file ('progress_status.txt' in the specified directory
 *  and writes "n% completed" etc in the file when  *  n% of a simulation has 
 *  been done, for integer n.
 */
class ProgressReporter
{
private:
    /*< Start time of the simulation */
    double mStartTime;
    /*< End time of the simulation */
    double mEndTime;
    /*< Progress status file */
    out_stream mpFile; 
    /*< Last percentage that was written */
    unsigned mLastPercentage;
       
public :

    /** 
     *  Constuctor saves times and opens output file ('progress_status.txt')
     */
    ProgressReporter(std::string outputDirectory, double startTime, double endTime)
        : mStartTime(startTime),
          mEndTime(endTime)
    {
        assert(startTime < endTime);
        
        // note we make sure we don't delete anything in the output directory
        OutputFileHandler handler(outputDirectory, false);
        mpFile = handler.OpenOutputFile("progress_status.txt");

        mLastPercentage = UINT_MAX;
    }
    
    ~ProgressReporter()
    {
        if(mLastPercentage!=100)
        {
            *mpFile << "100% completed" << std::endl;
        }
        *mpFile << "..done!" << std::endl;
        mpFile->close();
    }
    
    /**
     *  Calculates the percentage completed using the time given and the start and end
     *  time and prints to file if another percent has been done.
     */
    void Update(double currentTime)
    {
        unsigned percentage = floor( (currentTime - mStartTime)/(mEndTime - mStartTime)*100 );
        if(mLastPercentage==UINT_MAX || percentage > mLastPercentage)
        {
            *mpFile << percentage << "% completed" << std::endl;
            mLastPercentage = percentage;
        }
    }
    
    void PrintFinalising()
    {
        *mpFile << "Finalising.." << std::endl;
    }

    void PrintInitialising()
    {
        *mpFile << "Initialising.." << std::endl;
    }

};

#endif /*PROGRESSREPORTER_HPP_*/
