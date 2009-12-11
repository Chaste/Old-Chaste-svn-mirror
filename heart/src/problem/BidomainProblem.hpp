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

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <vector>

#include "AbstractCardiacProblem.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "Electrodes.hpp"
#include "BidomainPde.hpp"
#include "BidomainDg0Assembler.hpp"
#include "HeartRegionCodes.hpp"

/**
 * Class which specifies and solves a bidomain problem.
 *
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N),
 * where V_j is the voltage at node j and phi_j is the
 * extracellular potential at node j.
 */
template<unsigned DIM>
class BidomainProblem : public AbstractCardiacProblem<DIM,DIM, 2>
{

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save the member variables to an archive.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractCardiacProblem<DIM, DIM, 2> >(*this);
        archive & mpBidomainPde;
        //archive & mExtracelluarColumnId; // Created by InitialiseWriter, called from Solve
        archive & mRowForAverageOfPhiZeroed;
        archive & mHasBath;
        archive & mpElectrodes;
    }

    /**
     * Load the member variables from an archive.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacProblem<DIM, DIM, 2> >(*this);
        archive & mpBidomainPde;
        //archive & mExtracelluarColumnId; // Created by InitialiseWriter, called from Solve
        archive & mRowForAverageOfPhiZeroed;
        archive & mHasBath;
        archive & mpElectrodes;

        if (mHasBath)
        {
            // We only save element annotations, so annotate bath nodes from these
            AnalyseMeshForBath();
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    friend class TestBidomainWithBathAssembler;

protected:
    /** The bidomain PDE */
    BidomainPde<DIM>* mpBidomainPde;

    /** Nodes at which the extracellular voltage is fixed to zero (replicated) */
    std::vector<unsigned> mFixedExtracellularPotentialNodes;
   /** Used by the writer */
    unsigned mExtracelluarColumnId;
    /**
     * Another method of resolving the singularity in the bidomain equations.
     * Specifies a row of the matrix at which to impose an extra condition.
     */
    unsigned mRowForAverageOfPhiZeroed;

    /** Whether the mesh has a bath, ie whether this is a bath simulation */
    bool mHasBath;

    /** Electrodes used to provide a shock */
    Electrodes<DIM>* mpElectrodes;

    /**
     *  Create normal initial condition but overwrite V to zero for bath nodes, if
     *  there are any.
     */
    Vec CreateInitialCondition();

    /**
     * Annotate bath nodes with the correct region code, if a bath is present.
     * Will throw if #mHasBath is set but no bath is present in the mesh.
     */
    void AnalyseMeshForBath();

    /**
     *  We need to save the assembler that is being used to switch off the
     *  electrodes (by adding default boundary conditions to the assembler)
     */
    BidomainDg0Assembler<DIM,DIM>* mpAssembler;

    /** Create our bidomain PDE object */
    virtual AbstractCardiacPde<DIM> *CreateCardiacPde();

    /** Create a suitable bidomain assembler */
    virtual AbstractDynamicAssemblerMixin<DIM, DIM, 2>* CreateAssembler();

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     *   create cells.
     * @param hasBath Whether the simulation has a bath (if this is true, all elements with
     *   attribute = 1 will be set to be bath elements (the rest should have
     *   attribute = 0)).
     *
     */
    BidomainProblem(AbstractCardiacCellFactory<DIM>* pCellFactory, bool hasBath=false);

    /**
     * Constructor just used for archiving
     */
    BidomainProblem();

    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to
     *  zero. This does not necessarily have to be called. If it is not, phi_e
     *  is only defined up to a constant.
     *
     *  @param nodes  the nodes to be fixed.
     *
     *  @note currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes);

    /**
     * Set which row of the linear system should be used to enforce the
     * condition that the average of phi_e is zero.  If not called, this
     * condition will not be used.
     * 
     * @param node  the mesh node index giving the row at which to impose the constraint
     */
    void SetNodeForAverageOfPhiZeroed(unsigned node);

    /**
     *  Get the pde. Can only be called after Initialise()
     */
    BidomainPde<DIM>* GetBidomainPde();

    /**
     *  Print out time and max/min voltage/phi_e values at current time.
     * @param time  current time.
     */
    void WriteInfo(double time);

    /**
     * Define what variables are written to the primary results file.
     * Adds the extracellular potential.
     * @param extending  whether we are extending an existing results file
     */
    virtual void DefineWriterColumns(bool extending);

    /**
     * Write one timestep of output data to the primary results file.
     * Adds the extracellular potential to the results.
     * 
     * @param time  the current time
     * @param voltageVec  the solution vector to write
     */
    virtual void WriteOneStep(double time, Vec voltageVec);

    /**
     * Performs a series of checks before solving.
     * Checks that a suitable method of resolving the singularity is being used.
     */
    void PreSolveChecks();

    /**
     *  Set an electrode object (which provides boundary conditions). Only
     *  valid if there is a bath.
     * 
     * @param rElectrodes
     */
    void SetElectrodes(Electrodes<DIM>& rElectrodes);

    /**
     *  Called at end of each time step in the main time-loop in
     *  AbstractCardiacProblem::Solve(). Overloaded here to switch off
     *  the electrodes (if there are any).
     * 
     * @param time  the current time
     */
    void OnEndOfTimestep(double time);
};

#include "TemplatedExport.hpp" // Must be last
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainProblem);


#endif /*BIDOMAINPROBLEM_HPP_*/
