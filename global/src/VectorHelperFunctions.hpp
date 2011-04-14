/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef VECTORHELPERFUNCTIONS_HPP_
#define VECTORHELPERFUNCTIONS_HPP_

/**
 * @file
 *
 * A selection of helper functions to be able to access std::vector<double>
 * and CVODE's N_Vector types using the same interface.  These are used by
 * AbstractParameterisedSystem and some tests.
 */

#include <cassert>
#include <vector>

#ifdef CHASTE_CVODE
// CVODE headers
#include <nvector/nvector_serial.h>
#endif

/**
 * Helper function to get a vector component.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector to access
 * @param index  the index of the component to get
 */
template<typename VECTOR>
inline double GetVectorComponent(const VECTOR& rVec, unsigned index);

/**
 * Helper function to set a vector component.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector to modify
 * @param index  the index of the component to set
 * @param value  the new value
 */
template<typename VECTOR>
inline void SetVectorComponent(VECTOR& rVec, unsigned index, double value);

/**
 * Helper function to determine a vector's size.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector
 * @return  its size
 */
template<typename VECTOR>
inline unsigned GetVectorSize(const VECTOR& rVec);

/**
 * Helper function to initialise a vector to be empty/unset.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector
 */
template<typename VECTOR>
inline void InitialiseEmptyVector(VECTOR& rVec);

/**
 * Helper function to delete a vector.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector
 */
template<typename VECTOR>
inline void DeleteVector(VECTOR& rVec);

/**
 * Helper function to create a fresh copy of a vector.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector to copy
 */
template<typename VECTOR>
inline VECTOR CopyVector(VECTOR& rVec);


#ifdef CHASTE_CVODE

/**
 * A helper function to copy an N_Vector into a std::vector<double>.
 * Only exists if CHASTE_CVODE is defined.
 *
 * @param src  source vector
 * @param rDest  destination vector; will be resized and filled
 */
inline void CopyToStdVector(N_Vector src, std::vector<realtype>& rDest)
{
    // Check for no-op
    realtype* p_src = NV_DATA_S(src);
    if (p_src == &(rDest[0])) return;
    // Set dest size
    long size = NV_LENGTH_S(src);
    rDest.resize(size);
    // Copy data
    for (long i=0; i<size; i++)
    {
        rDest[i] = p_src[i];
    }
}


/**
 * A helper function to copy a std::vector<double> into an N_Vector.
 * Only exists if CHASTE_CVODE is defined.
 *
 * @param rSrc  source vector
 * @param dest  destination vector; must exist and be the correct size
 */
inline void CopyFromStdVector(const std::vector<realtype>& rSrc, N_Vector dest)
{
    // Check for no-op
    realtype* p_dest = NV_DATA_S(dest);
    if (p_dest == &(rSrc[0])) return;
    // Check dest size
    long size = NV_LENGTH_S(dest);
    assert(size == (long)rSrc.size());
    // Copy data
    for (long i=0; i<size; i++)
    {
        p_dest[i] = rSrc[i];
    }
}

#endif // CHASTE_CVODE


//
// Specialisations for std::vector<double>
//

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @param index
 */
template<>
inline double GetVectorComponent(const std::vector<double>& rVec, unsigned index)
{
    assert(index < rVec.size());
    return rVec[index];
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @param index
 * @param value
 */
template<>
inline void SetVectorComponent(std::vector<double>& rVec, unsigned index, double value)
{
    assert(index < rVec.size());
    rVec[index] = value;
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 */
template<>
inline unsigned GetVectorSize(const std::vector<double>& rVec)
{
    return rVec.size();
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 */
template<>
inline void InitialiseEmptyVector(std::vector<double>& rVec)
{
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 */
template<>
inline void DeleteVector(std::vector<double>& rVec)
{
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 */
template<>
inline std::vector<double> CopyVector(std::vector<double>& rVec)
{
    return rVec;
}

//
// Specialisations for N_Vector
//

#ifdef CHASTE_CVODE

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @param index
 */
template<>
inline double GetVectorComponent(const N_Vector& rVec, unsigned index)
{
    assert(rVec != NULL);
    return NV_Ith_S(rVec, index);
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @param index
 * @param value
 */
template<>
inline void SetVectorComponent(N_Vector& rVec, unsigned index, double value)
{
    assert(rVec != NULL);
    NV_Ith_S(rVec, index) = value;
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 */
template<>
inline unsigned GetVectorSize(const N_Vector& rVec)
{
    assert(rVec != NULL);
    return NV_LENGTH_S(rVec);
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 */
template<>
inline void InitialiseEmptyVector(N_Vector& rVec)
{
    rVec = NULL;
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 */
template<>
inline void DeleteVector(N_Vector& rVec)
{
    if (rVec)
    {
        rVec->ops->nvdestroy(rVec);
        rVec = NULL;
    }
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 */
template<>
inline N_Vector CopyVector(N_Vector& rVec)
{
    N_Vector copy = NULL;
    if (rVec)
    {
        copy = N_VClone(rVec);
        unsigned size = NV_LENGTH_S(rVec);
        for (unsigned i=0; i<size; i++)
        {
            NV_Ith_S(copy, i) = NV_Ith_S(rVec, i);
        }
    }
    return copy;
}
#endif // CHASTE_CVODE

//
// End of helper functions
//


#endif /*VECTORHELPERFUNCTIONS_HPP_*/
