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
 * See tutorials for usage.
 */
template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractCardiacProblem
{
friend class TestBidomainWithBathAssembler;
    
protected:
    /** Meshes can be read from file or instantiated and passed directly to this 
     *  class, this is for the former */ 
    std::string mMeshFilename;

    /** If this is set, the nodes for each processor are read */
    std::string mNodesPerProcessorFilename;

    /** Data is not written if output directory or output file prefix are not set*/
    std::string  mOutputDirectory, mOutputFilenamePrefix;

    /**
     *  Whether to use matrix-based assembly of the RHS vector (much more efficient). 
     *  True by default
     */ 
    bool mUseMatrixBasedRhsAssembly;
    /** Whether this problem class has created the mesh itself, as opposed to being given it */
    bool mAllocatedMemoryForMesh;
    /** Whether to print some statistics (max/min voltage) to screen during the simulation */
    bool mWriteInfo;
    /** Whether to write any output at all */
    bool mPrintOutput;
    /** Whether to convert the output from HDF5 to meshalyzer readable format */ 
    bool mCallChaste2Meshalyzer;

    /** If only outputing voltage for selected nodes, which nodes to output at */
    std::vector<unsigned> mNodesToOutput;

    /** Used by the writer */
    unsigned mVoltageColumnId;
    /** Used by the writer */
    unsigned mTimeColumnId;
    /** Used by the writer */
    unsigned mNodeColumnId;

    /** The monodomain or bidomain pde */
    AbstractCardiacPde<SPACE_DIM>* mpCardiacPde;

    /** Boundary conditions container used in the simulation */
    BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditionsContainer;
    /** It is convenient to also have a separate variable for default (zero-Neumann) boundary conditions */
    BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpDefaultBoundaryConditionsContainer;
    /** The PDE solver */
    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpAssembler;
    /** The cell factory creates the cells for each node */
    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    /** The mesh. Can either by passed in, or the mesh filename can be set */
    AbstractMesh<SPACE_DIM,SPACE_DIM>* mpMesh;

    /** The current solution vector, of the form [V_0 .. V_N ] for monodomain and 
     *  [V_0 phi_0 .. V_N phi_N] for bidomain */
    Vec mSolution; 

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
    // ///\todo CardiacElectroMechanicsProblem should be a friend, but not sure
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

    /** 
     *  Set a file from which the nodes for each processor are read
     */ 
    void SetNodesPerProcessorFilename(const std::string& filename);
    
    /**
     *  Set the boundary conditions container.
     *  @param *pbcc is a pointer to a boundary conditions container
     */
    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM> *bcc);
    
    /**
     *  Performs a series of checks before solving.
     *  It checks whether the cardiac pde has been defined, 
     *  whether the simulation time is greater than zero and 
     *  whether the output directory is specified (or the output is set not to be produced).
     *  It throws exceptions if any of the above checks fails.
     */
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

    /** This only needs to be called if a mesh filename has not been set */ 
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
    
    /**
     *  First performs some checks by calling  the PreSolveChecks method. 
     *  It creates an assembler to which it passes the boundary conditions specified by the user
     *  (otherwise it passes the defauls bcc). 
     *   It then calls the Solve method in the assembler class.
     *   It also handles the output, if necessary.
     */
    void Solve();
    
    /**
     * Closes the files where the solution is stored and,
     * if specified so (as it is by default), converts the output to Meshalyzer format 
     * by calling the WriteFilesUsingMesh method in the MeshalyzerWriter class.
     */
    void CloseFilesAndPostProcess();


    virtual void WriteInfo(double time)=0;

    virtual void DefineWriterColumns();

    virtual void WriteOneStep(double time, Vec voltageVec);
    
    /**
     * It creates and initialises the hdf writer from the Hdf5DataWriter class. 
     * It passes the output directory (mOutputDirectory) and file name (mOutputFilenamePrefix) to it.
     * It is called by Solve(), if the output needs to be generated.
     */
    void InitialiseWriter();
    
    /**
     * Specifies which nodes in the mesh to output.
     * @param &nodesToOutput is a reference to a vector with the indexes of the nodes
     * where the output is desired. 
     * If empty, the output will be for all the nodes in the mesh.
     */
    void SetOutputNodes(std::vector<unsigned> &nodesToOutput);
    
    Hdf5DataReader GetDataReader();
    
    /** 
     *  Whether to use matrix-based RHS assembly or not
     */    
    void UseMatrixBasedRhsAssembly(bool usematrix=true);
    
    /**
     *  Called at end of each time step in the main time-loop in 
     *  Solve(). Empty implementation but can be overloaded by child 
     *  classes
     */
    virtual void OnEndOfTimestep(double time)
    {}
    
};
#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
