/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef PARALLELCOLUMNDATAWRITER_HPP_
#define PARALLELCOLUMNDATAWRITER_HPP_


#include "ColumnDataWriter.hpp"
#include "DistributedVector.hpp"

#include <petscvec.h>

/**
 * A parallelised column data writer class.
 */
class ParallelColumnDataWriter  : public ColumnDataWriter
{
private:

    bool mIsParallel;     /**< Set to true in constructor if running in parallel*/
    bool mAmMaster;       /**< Set to true in constructor for process is the rank 0 process*/
    Vec mConcentrated;    /**< Vector to hold concentrated copy of distributed vector on the master process*/
    VecScatter mToMaster; /**< variable holding information for concentrating a vector*/

public:

    /**
     * Constructor.
     * 
     * @param directory  the directory in which to write the data to file
     * @param baseName  the name of the file in which to write the data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     */
    ParallelColumnDataWriter(std::string directory, std::string baseName, bool cleanDirectory=true);

    /**
     * Destructor.
     */
    virtual ~ParallelColumnDataWriter();

    bool AmMaster() const;
    void PutVector(int variableID, Vec PetscVector);
    void PutVectorStripe(int variableId, DistributedVector::Stripe stripe);

    /**
     * Input the variable value to the output file or ancillary file
     * 
     * @param variableID
     * @paramvariableValue
     * @param dimensionPosition  The position in column (defaults to -1). This is required if 
     *      there is a fixed dimension, and will be the position along that dimension
     */
    void PutVariable(int variableID, double variableValue, long dimensionPosition = -1);

    /**
     * End the define mode of the DataWriter.
     */
    void EndDefineMode();

    /**
     * Advance along the unlimited dimension. Normally this will be called
     * when all variables in a row have been input.
     */
    void AdvanceAlongUnlimitedDimension();

    /**
     * Close any open files.
     */
    void Close();
};

#endif /*PARALLELCOLUMNDATAWRITER_HPP_*/
