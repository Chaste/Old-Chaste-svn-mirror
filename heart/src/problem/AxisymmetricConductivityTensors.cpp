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

#include <vector>
#include "UblasIncludes.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "Exception.hpp"

template<unsigned SPACE_DIM>
void AxisymmetricConductivityTensors<SPACE_DIM>::ReadOrientationVectorFromFile (c_vector<double,SPACE_DIM>& rOrientVector)
{
    std::vector<double> tokens;

    unsigned num_elems = this->GetTokensAtNextLine(tokens);

    if (num_elems != 3u)
    {
        this->CloseFibreOrientationFile();
        EXCEPTION("Axisymmetric media defined. Fibre orientation file should contain 3 values per element");
    }

    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        rOrientVector[i] = tokens[i];
    }
}

template<unsigned SPACE_DIM>
AxisymmetricConductivityTensors<SPACE_DIM>::AxisymmetricConductivityTensors()
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("Axisymmetric anisotropic conductivity only makes sense in 3D");
    }
}

template<unsigned SPACE_DIM>
void AxisymmetricConductivityTensors<SPACE_DIM>::SetConstantConductivities(c_vector<double, 3> constantConductivities)
{
    //assert(SPACE_DIM == 3);//Otherwise constructor would have thrown
    if (constantConductivities[1] != constantConductivities[2])
    {
        EXCEPTION("Axisymmetric media defined: transversal and normal conductivities should have the same value");
    }

    this->mUseNonConstantConductivities = false;
    this->mConstantConductivities = constantConductivities;
}

template<unsigned SPACE_DIM>
void AxisymmetricConductivityTensors<SPACE_DIM>::Init() throw (Exception)
{
    if (!this->mUseNonConstantConductivities && !this->mUseFibreOrientation)
    {
        // Constant tensor for every element
        c_matrix<double, SPACE_DIM, SPACE_DIM> conductivity_matrix(zero_matrix<double>(SPACE_DIM,SPACE_DIM));

        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            assert(this->mConstantConductivities(dim) != DBL_MAX);
            conductivity_matrix(dim,dim) = this->mConstantConductivities(dim);
        }
        this->mTensors.push_back(conductivity_matrix);
    }
    else
    {
        c_vector<double,SPACE_DIM> fibre_vector((zero_vector<double>(SPACE_DIM)));
        fibre_vector[0]=1.0;

        if (this->mUseFibreOrientation)
        {
            this->OpenFibreOrientationFile();
            this->mNumElements = this->GetNumElementsFromFile();
        }
        else
        {
            this->mNumElements = this->mpNonConstantConductivities->size();
        }

        // reserve() allocates all the memory at once, more efficient than relying
        // on the automatic reallocation scheme.
        this->mTensors.reserve(this->mNumElements);

        c_matrix<double, SPACE_DIM, SPACE_DIM> conductivity_matrix(zero_matrix<double>(SPACE_DIM,SPACE_DIM));

        for (unsigned element_index=0; element_index<this->mNumElements; element_index++)
        {
            /*
             *  For every element of the mesh we compute its tensor like (from
             * "Laminar Arrangement of VentricularMyocites Influences Electrical
             * Behavior of the Heart", Darren et al. 2007):
             *
             *                         [g_f  0   0 ] [a_f']
             *  tensor = [a_f a_l a_n] [ 0  g_l  0 ] [a_l']
             *                         [ 0   0  g_n] [a_n']
             *
             *              [x_i]
             *  where a_i = [y_i], i={f,l,n} are read from the fibre orientation file and
             *              [z_i]
             *
             *  g_f = fibre/longitudinal conductivity (constant or element specific)
             *  g_l = laminar/transverse conductivity (constant or element specific)
             *  g_n = normal conductivity (constant or element specific)
             *
             *
             *  For axisymmetric anisotropic media (g_l = g_n) we can simplify previous expression to
             *
             *
             *  tensor = g_l * I + (g_f - g_l) * a_f * a_f'
             *
             */

            if (this->mUseNonConstantConductivities)
            {
                for (unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    conductivity_matrix(dim,dim) = (*this->mpNonConstantConductivities)[element_index][dim];
                }
            }
            else
            {
                for (unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    assert(this->mConstantConductivities(dim) != DBL_MAX);
                    conductivity_matrix(dim,dim) = this->mConstantConductivities(dim);
                }
            }


            if (this->mUseFibreOrientation)
            {
                ReadOrientationVectorFromFile(fibre_vector);
            }

            this->mTensors.push_back( conductivity_matrix(1,1) * identity_matrix<double>(SPACE_DIM) +
                                      (conductivity_matrix(0,0) - conductivity_matrix(1,1)) * outer_prod(fibre_vector,fibre_vector));
        }

        if (this->mUseFibreOrientation)
        {
            this->CloseFibreOrientationFile();
        }
    }

    this->mInitialised = true;
}



/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

// only makes sense in 3d, but we need the other to compile
// AbstractCardiacPde and BidomainTissue.
template class AxisymmetricConductivityTensors<1>;
template class AxisymmetricConductivityTensors<2>;
template class AxisymmetricConductivityTensors<3>;
