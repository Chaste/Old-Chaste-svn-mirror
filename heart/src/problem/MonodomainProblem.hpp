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


#ifndef MONODOMAINPROBLEM_HPP_
#define MONODOMAINPROBLEM_HPP_

#include "AbstractCardiacProblem.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "MonodomainPde.hpp"


/**
 * Class which specifies and solves a monodomain problem.
 */
template<unsigned DIM>
class MonodomainProblem : public AbstractCardiacProblem<DIM, 1>
{
protected:
    MonodomainPde<DIM>* mpMonodomainPde;

public:
    AbstractCardiacPde<DIM>* CreateCardiacPde();

    AbstractDynamicAssemblerMixin<DIM, DIM, 1>* CreateAssembler();

public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    MonodomainProblem(AbstractCardiacCellFactory<DIM>* pCellFactory);

    /**
     * Destructor
     */
    ~MonodomainProblem();

    MonodomainPde<DIM> * GetMonodomainPde();

    /**
     *  Print out time and max/min voltage values at current time.
     */
    void WriteInfo(double time);
};

#endif /*MONODOMAINPROBLEM_HPP_*/
