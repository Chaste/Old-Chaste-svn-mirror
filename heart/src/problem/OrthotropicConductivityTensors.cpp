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

#include "OrthotropicConductivityTensors.hpp"
#include "Exception.hpp"

template<unsigned SPACE_DIM>
void OrthotropicConductivityTensors<SPACE_DIM>::ReadOrientationMatrixFromFile (c_matrix<double,SPACE_DIM,SPACE_DIM>& orientMatrix)
{
    std::vector<double> items;
    unsigned num_elems = this->GetTokensAtNextLine(items);

    if (num_elems != SPACE_DIM*SPACE_DIM)
    {
        this->CloseFibreOrientationFile();
        EXCEPTION("Orthotropic media defined. Number of vectors in fibre orientation file and size of them should match SPACE_DIM");
    }

    for (unsigned vector_index=0; vector_index<SPACE_DIM; vector_index++)
    {
        for (unsigned component_index=0; component_index<SPACE_DIM; component_index++)
        {
            orientMatrix(component_index, vector_index) = items[vector_index*SPACE_DIM + component_index];
        }
    }
}

template<unsigned SPACE_DIM>
void OrthotropicConductivityTensors<SPACE_DIM>::Init() throw (Exception)
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
        c_matrix<double,SPACE_DIM,SPACE_DIM> orientation_matrix((identity_matrix<double>(SPACE_DIM)));

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
                ReadOrientationMatrixFromFile(orientation_matrix);
            }

            c_matrix<double,SPACE_DIM,SPACE_DIM> temp;
            noalias(temp) = prod(orientation_matrix, conductivity_matrix);
            this->mTensors.push_back( prod(temp, trans(orientation_matrix) ) );
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

template class OrthotropicConductivityTensors<1>;
template class OrthotropicConductivityTensors<2>;
template class OrthotropicConductivityTensors<3>;
