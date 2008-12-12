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
#ifndef CELLWISENUTRIENTSINKPDE_HPP_
#define CELLWISENUTRIENTSINKPDE_HPP_

#include "MeshBasedTissue.hpp"
#include "AbstractLinearEllipticPde.hpp"

/**
 *  A nutrient PDE which has a sink at each non-apoptotic cell.
 */
template<unsigned DIM>
class CellwiseNutrientSinkPde : public AbstractLinearEllipticPde<DIM>
{
private:

    MeshBasedTissue<DIM>& mrTissue;

    double mCoefficient;

public:

    CellwiseNutrientSinkPde(MeshBasedTissue<DIM>& rTissue, double coefficient);

    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x);

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x, Element<DIM,DIM>*);

    double ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode);

    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& );

};

#endif /*CELLWISENUTRIENTSINKPDE_HPP_*/
