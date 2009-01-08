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


#ifndef PETSCTOOLS_HPP_
#define PETSCTOOLS_HPP_


#include <vector>
#include <iostream>
#include <cassert>
#include <cstring> //For strcmp etc. Needed in gcc-4.3

#ifndef SPECIAL_SERIAL
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"
#endif //SPECIAL_SERIAL

#define EXIT_IF_PARALLEL if(!PetscTools::IsSequential()){TS_TRACE("This test does not pass in parallel yet.");return;}
#define EXIT_IF_SEQUENTIAL if(PetscTools::IsSequential()){TS_TRACE("This test is not meant to be executed in sequential.");return;}

/**
 * A helper class of static methods.
 */
class PetscTools
{
private:
    /** Whether PETSc has been initialised. */
    static bool mPetscIsInitialised;
    
    /** The total number of processors. */
    static unsigned mNumProcessors;
    
    /** Which processors we are. */
    static unsigned mRank;
    
public:
    
    /**
     * Reset our cached values: whether PETSc is initialised,
     * how many processors there are, and which one we are.
     */
    static void ResetCache();
    
    /**
     * Just returns whether there is one process or not.
     */
    static bool IsSequential();
    
    /**
     *  Returns total number of processors
     */
    static unsigned NumProcs();
    
    /**
     * Return our rank.
     * 
     * If PETSc has not been initialized, returns 0.
     */
    static unsigned GetMyRank();
    
    /**
     * Just returns whether it is the master process or not.
     * 
     * If not running in parallel, always returns true.
     */
    static bool AmMaster();
    
    /**
     * If MPI is set up, perform a barrier synchronisation.
     * If not, it's a noop.
     */
    static void Barrier();

#ifndef SPECIAL_SERIAL
    /**
     *  Create a vector of the specified size. SetFromOptions is called.
     */
    static Vec CreateVec(int size);

    /**
     *  Create a vector of the specified size with all values set to be the given
     *  constant. SetFromOptions is called.
     */
    static Vec CreateVec(int size, double value);

    /**
     *  Create a Vec from the given data.
     */
    static Vec CreateVec(std::vector<double> data);

    /**
     * Set up a matrix - set the size using the given parameters, the type (default MATMPIAIJ). The
     * number of local rows and columns is by default PETSC_DECIDE. SetFromOptions is called.
     *
     * @param maxColsPerRow The maximum number of non zeros per row. This value is problem dependent.
     *     An upper bound is (3^ELEMENT_DIM) * PROBLEM_DIM. The default value (3D bidomain problem)
     *     should be big enough for any of the problems being solved.
     */
    static void SetupMat(Mat& rMat, int numRows, int numColumns,
                         MatType matType=(MatType) MATMPIAIJ,
                         int numLocalRows=PETSC_DECIDE,
                         int numLocalColumns=PETSC_DECIDE,
                         int maxColsPerRow=54);

    /**
     * Ensure exceptions are handled cleanly in parallel code, by causing all processes to
     * throw if any one does.
     *
     * @param flag is set to true if this process has thrown.
     */
    static void ReplicateException(bool flag);

    /**
     *  Another helper method to get a single value from a vector
     *  in 1 line than Petsc's usual 4 or 5. DOES NOT check that
     *  the requested component is local, DOES do bound-checking.
     * 
     * \todo Think if there is an efficient compromise between this method and
     * the full weight of ReplicatableVector (broadcast single values to all processors).
     *  How do you know who has the value?
     * 
     * 
     */
//    static double GetVecValue(Vec vec, unsigned index)
//    {
//        assert(vec);
//        PetscInt size;
//        VecGetSize(vec, &size);
//        assert((int)index<size);
//
//        double* p_data;
//        VecGetArray(vec, &p_data);
//        double ret = p_data[(int)index];
//        VecRestoreArray(vec, &p_data);
//
//        return ret;
//    }

#endif //SPECIAL_SERIAL

};



#endif /*PETSCTOOLS_HPP_*/
