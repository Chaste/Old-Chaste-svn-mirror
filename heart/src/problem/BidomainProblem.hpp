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


#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_

#include <vector>

#include "AbstractCardiacProblem.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "Electrodes.hpp"
#include "BidomainPde.hpp"

/**
 * Class which specifies and solves a bidomain problem.
 *
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N),
 * where V_j is the voltage at node j and phi_j is the
 * extracellular potential at node j.
 */
template<unsigned SPACE_DIM>
class BidomainProblem : public AbstractCardiacProblem<SPACE_DIM, 2>
{

friend class TestBidomainWithBathAssembler;    
    
protected:
    BidomainPde<SPACE_DIM>* mpBidomainPde;

private:
    std::vector<unsigned> mFixedExtracellularPotentialNodes; /** nodes at which the extracellular voltage is fixed to zero (replicated) */
    unsigned mExtracelluarColumnId;
    unsigned mRowMeanPhiEZero;
    bool mHasBath;

    /**
     *  Create normal initial condition but overwrite V to zero for bath nodes, if 
     *  there are any.
     */
    Vec CreateInitialCondition();
    
    /**
     * Annotate bath nodes with the correct region code, if a bath is present.
     * Will throw if mHasBath is set but no bath is present in the mesh.
     */
    void AnalyseMeshForBath();

protected:
    AbstractCardiacPde<SPACE_DIM> *CreateCardiacPde();

    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, 2>* CreateAssembler();

public:
    /**
     * Constructor 
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    BidomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, bool hasBath=false);

    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to
     *  zero. This does not necessarily have to be called. If it is not, phi_e
     *  is only defined up to a constant.
     *
     *  @param the nodes to be fixed.
     *
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes);

    /**
     * Set which row of the linear system should be used to enforce the
     * condition that the mean of phi_e is zero.  If not called, this
     * condition will not be used.
     */
    void SetRowForMeanPhiEToZero(unsigned rowMeanPhiEZero);
    
    /**
     *  Get the pde. Can only be called after Initialise()
     */
    BidomainPde<SPACE_DIM>* GetBidomainPde();

    /**
     *  Print out time and max/min voltage/phi_e values at current time.
     */
    void WriteInfo(double time);

    virtual void DefineWriterColumns();

    virtual void WriteOneStep(double time, Vec voltageVec);
    
    void PreSolveChecks();
    
    void SetElectrodes(Electrodes<SPACE_DIM>& rElectrodes);
};


#endif /*BIDOMAINPROBLEM_HPP_*/
