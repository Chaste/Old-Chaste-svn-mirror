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

#ifndef DISTRIBUTEDVECTORFACTORY_HPP_
#define DISTRIBUTEDVECTORFACTORY_HPP_

#include <boost/serialization/access.hpp>
#include <petscvec.h>
#include <cassert>
#include "DistributedVector.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * Factory for creating PETSc vectors distributed across processes.
 *
 * Replacement for the vector creation portions of DistributedVector (which
 * was implemented using static methods and data), the factory class allows
 * several patterns of PETSc vector length (and distributions among
 * processes) to co-exist.
 *
 * All vectors created by a factory instance will have the same base size
 * and parallelisation pattern.
 */
 class DistributedVectorFactory
{
private:
    // Data global to all vectors created by this factory.
    /** The first entry owned by the current processor. */
    unsigned mLo;
    /** One above the last entry owned by the current processor. */
    unsigned mHi;
    /** The problem size, i.e. the number of nodes in the mesh (the number of unknowns may be larger in a Stripe). */
    unsigned mProblemSize;
    /** Whether we've checked that PETSc is initialised. */
    bool mPetscStatusKnown;
    
    /**
     * Double check (in debug code) that PETSc has been initialised properly
     */
    void CheckForPetsc();  
    
    /**
     * Helper method for the constructors
     * 
     * @param vec the sample PETSc vector from which to calculate ownerships 
     */
    void CalculateOwnership(Vec vec);
    
    /** Needed for serialization. */
    friend class boost::serialization::access;
    
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Nothing to do - all done in load_construct_data
    }    
    
    
public:
    /**
     * Set the problem with an existing PETSc vector -- must have stride=1.
     *
     * @param vec is a PETSc vector which we want to use as the pattern for future vectors produced by this factory
     */
    DistributedVectorFactory(Vec vec);

    /**
     * Set the problem size specifying distribution over local processor.
     *
     * @param size
     * @param local - default to PETSc's default
     */
    DistributedVectorFactory(unsigned size, PetscInt local=PETSC_DECIDE);


    /**
     * Create a PETSc vector of the problem size
     * 
     * @return the created vector
     */    
    Vec CreateVec();
    
    /**
     * Create a striped PETSc vector of size: stride * problem size
     *
     * @param stride
     */
    Vec CreateVec(unsigned stride);
    
    
    /**
     * Create a distributed vector which wraps a given petsc vector
     * 
     * @param vec is the vector
     * @return the distributed vector 
     */
    DistributedVector CreateDistributedVector(Vec vec);
    
    /**
     * Create a distributed vector stripe which wraps a given petsc vector
     * 
     * @param vec is the vector
     * @return the distributed vector 
     */
//    DistributedVector::Stripe CreateStripe(Vec vec, unsigned stripe);
    
    /**
     * @return The number of elements in the vector owned by the local process
     */
    unsigned GetLocalOwnership() const
    {
        return mHi - mLo;
    }
    
    /**
     * @return mHi - The next index above the top one owned by the process.
     */
    unsigned GetHigh() const
    {
        return mHi;
    }
    
    /**
     * @return mLo - The lowest index owned by the process.
     */
    unsigned GetLow() const
    {
        return mLo;
    }

    /**
     * @return The number of elements in the vector
     */
    unsigned GetSize() const
    {
        return mProblemSize;
    }
    
};
// Declare identifier for the serializer
BOOST_CLASS_EXPORT(DistributedVectorFactory);

namespace boost
{
namespace serialization
{

template<class Archive>
inline void save_construct_data(
    Archive & ar, const DistributedVectorFactory * t, const unsigned int file_version)
{
    unsigned num_procs, lo, hi, size;
    hi = t->GetHigh();
    ar << hi;
    lo = t->GetLow();
    ar << lo;
    size = t->GetSize();
    ar << size;
    num_procs = PetscTools::GetNumProcs();
    ar << num_procs;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, DistributedVectorFactory * t, const unsigned int file_version)
{
    unsigned num_procs, lo, hi, size;
    ar >> hi;
    ar >> lo;
    ar >> size;
    ar >> num_procs;
    if (num_procs != PetscTools::GetNumProcs())
    {
        EXCEPTION("This archive was written for a different number of processors");
    }
    ::new(t)DistributedVectorFactory(size, hi-lo);
}

}
} // namespace ...

#endif /*DISTRIBUTEDVECTORFACTORY_HPP_*/
