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


#ifndef PETSCTOOLS_HPP_
#define PETSCTOOLS_HPP_

#include <string>
#include <vector>

#ifndef SPECIAL_SERIAL
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
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
    
    /** Private method makes sure that (if this is the first use within a test) then PETSc has been probed */
    static inline void CheckCache()
    {
        if (mNumProcessors == 0)
        {
            ResetCache();
        }
    }

public:

    /**
     * As a convention, we consider processor 0 the master process
     */
    static const unsigned MASTER_RANK=0;

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
    static unsigned GetNumProcs();

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
     * Just returns whether it is the right-most process or not.
     *
     * If not running in parallel, always returns true.
     */
    static bool AmTopMost();

    /**
     * If MPI is set up, perform a barrier synchronisation.
     * If not, it's a noop.
     */
    static void Barrier();

#ifndef SPECIAL_SERIAL
    /**
     * Create a vector of the specified size. SetFromOptions is called.
     *
     * @param size  the size of the vector
     */
    static Vec CreateVec(int size);

    /**
     * Create a vector of the specified size with all values set to be the given
     * constant. SetFromOptions is called.
     *
     * @param size  the size of the vector
     * @param value  the value to set each entry
     */
    static Vec CreateVec(int size, double value);

    /**
     * Create a Vec from the given data.
     *
     * @param data  some data
     */
    static Vec CreateVec(std::vector<double> data);

    /**
     * Set up a matrix - set the size using the given parameters, the type (default MATMPIAIJ). The
     * number of local rows and columns is by default PETSC_DECIDE. SetFromOptions is called.
     *
     * @param rMat the matrix
     * @param numRows the number of rows in the matrix
     * @param numColumns the number of columns in the matrix
     * @param matType the matrix type (defaults to MATMPIAIJ)
     * @param numLocalRows the number of local rows (detaults to PETSC_DECIDE)
     * @param numLocalColumns the number of local columns (detaults to PETSC_DECIDE)
     * @param maxColsPerRowIfMatMpiAij The maximum number of non zeros per row. This value is problem dependent.
     *   Since the call to set this depends on the matrix-type (eg MatMPIAIJSetPreallocation/MatSeqAIJSetPreallocation),
     *   preallocation using this value is done only if the matrix-type is MATMPIAIJ (the default). WITH OTHER
     *   TYPES OF MATRIX NO PREALLOCATION IS DONE AND YOU MUST PREALLOCATE MANUALLY (by calling the appropriate method)!
     */
    static void SetupMat(Mat& rMat, int numRows, int numColumns,
                         MatType matType=(MatType) MATMPIAIJ,
                         int numLocalRows=PETSC_DECIDE,
                         int numLocalColumns=PETSC_DECIDE,
                         int maxColsPerRowIfMatMpiAij=54 /* see doxygen comment! */);

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

    /**
     * Dumps a given Petsc object to disk.
     *
     * @param rMat a matrix
     * @param rOutputFileFullPath where to dump the matrix to disk
     */
    static void DumpPetscObject(const Mat& rMat, const std::string& rOutputFileFullPath);

    /**
     * Dumps a given Petsc object to disk.
     *
     * @param rVec a vector
     * @param rOutputFileFullPath where to dump the vector to disk
     */
    static void DumpPetscObject(const Vec& rVec, const std::string& rOutputFileFullPath);

    /**
     * Read a previously dumped Petsc object from disk.
     *
     * @param rMat a matrix
     * @param rOutputFileFullPath where to read the matrix from
     */
    static void ReadPetscObject(Mat& rMat, const std::string& rOutputFileFullPath);

    /**
     * Read a previously dumped Petsc object from disk.
     *
     * @param rVec a vector
     * @param rOutputFileFullPath where to read the matrix from
     */
    static void ReadPetscObject(Vec& rVec, const std::string& rOutputFileFullPath);

#endif //SPECIAL_SERIAL

};

#endif /*PETSCTOOLS_HPP_*/
