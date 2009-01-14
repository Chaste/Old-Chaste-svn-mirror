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


#ifndef ABSTRACTCARDIACPROBLEM_HPP_
#define ABSTRACTCARDIACPROBLEM_HPP_

#include <string>
#include <vector>

#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractMesh.hpp"

#include "BoundaryConditionsContainer.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5DataWriter.hpp"

/**
 * Base class for cardiac problems; contains code generic to both mono- and bidomain.
 * 
 * \todo further documentation, including of member variables.
 */
template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractCardiacProblem
{
friend class TestBidomainWithBathAssembler;
    
protected:
    std::string mMeshFilename;
    std::string mNodesPerProcessorFilename;

    /** data is not written if output directory or output file prefix are not set*/
    std::string  mOutputDirectory, mOutputFilenamePrefix;

    bool mUseMatrixBasedRhsAssembly;
    bool mAllocatedMemoryForMesh;
    bool mWriteInfo;
    bool mPrintOutput;
    bool mCallChaste2Meshalyzer;

    std::vector<unsigned> mNodesToOutput;

    unsigned mVoltageColumnId;
    unsigned mTimeColumnId;
    unsigned mNodeColumnId;

    AbstractCardiacPde<SPACE_DIM>* mpCardiacPde;

    BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditionsContainer;
    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpAssembler;

    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    AbstractMesh<SPACE_DIM,SPACE_DIM>* mpMesh;

    Vec mSolution; // Current solution

    /**
     * Subclasses must override this method to create a PDE object of the appropriate type.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractCardiacPde<SPACE_DIM>* CreateCardiacPde() =0;

    /**
     * Subclasses must override this method to create a suitable assembler object.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* CreateAssembler() =0;

public:
    // This (and things in MonodomainProblem) being public are hacks for
    // CardiacElectroMechanicsProblem to work.
    // TODO CardiacElectroMechanicsProblem should be a friend, but not sure
    // how to get friends to work when both friends are templated and abstract.
    Hdf5DataWriter* mpWriter;

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    AbstractCardiacProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory);

    virtual ~AbstractCardiacProblem();

    /*
     *  Initialise the system. Must be called before Solve()
     */
    void Initialise();

    void SetNodesPerProcessorFilename(const std::string& filename);

    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM> *bcc);

    virtual void PreSolveChecks();

    /**
     * Perhaps this should be a method of AbstractCardiacPde??
     * This is virtual so BidomainProblem can overwrite V to zero for bath nodes, if 
     * there are any.
     */
    virtual Vec CreateInitialCondition();

    /** 
     *  Set whether to call the Chaste2Meshalyzer script.
     *  This script gets everything ready to visualize the results with meshalyser
     *  and is useful in testing. By default the script is called.
     *  In performance testing for example it desirable to disable the script.
     */
    void ConvertOutputToMeshalyzerFormat(bool call = true);

    void SetMesh(AbstractMesh<SPACE_DIM,SPACE_DIM>* pMesh);

    /**
     *  Set whether the simulation will generate results files.
     */
    void PrintOutput(bool rPrintOutput);

    /**
     *  Set whether extra info will be written to stdout during computation.
     */
    void SetWriteInfo(bool writeInfo = true);

    /**
     *  Get the final solution vector. This vector is distributed over all processes.
     *
     *  In case of Bidomain, this is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N).
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     *
     *  Use with caution since we don't want to alter the state of the PETSc vector.
     */
    Vec GetSolution();

    AbstractMesh<SPACE_DIM,SPACE_DIM> & rGetMesh();

    AbstractCardiacPde<SPACE_DIM>* GetPde();

    void Solve();

    void CloseFilesAndPostProcess();


    virtual void WriteInfo(double time) =0;

    virtual void DefineWriterColumns();

    virtual void WriteOneStep(double time, Vec voltageVec);

    void InitialiseWriter();

    void SetOutputNodes(std::vector<unsigned> &nodesToOutput);
    
    Hdf5DataReader GetDataReader();
    
    /** 
     *  Whether to use matrix-based RHS assembly or not
     */    
    void UseMatrixBasedRhsAssembly(bool usematrix=true);
    
};
#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
