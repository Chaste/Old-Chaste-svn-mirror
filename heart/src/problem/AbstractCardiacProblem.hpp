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
#include <cassert>
#include <climits> // Work around a boost bug - see #1024.
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "ClassIsAbstract.hpp"
#include "AbstractCardiacCell.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "DistributedVectorFactory.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5DataWriter.hpp"

/*
 * Archiving extravaganza: 
 * 
 * We archive mesh and pde through a pointer to an abstract class. All the potential concrete
 * classes need to be included here, so they are registered with boost. If not, boost won't be
 * able to find the archiving methods of the concrete class and will throw the following
 * exception:
 * 
 *       terminate called after throwing an instance of 'boost::archive::archive_exception'
 *       what():  unregistered class
 * 
 * No member variable is defined to be of any of these clases, removing them won't
 * produce any compiler error. The exception above will occur at runtime.
 * 
 * This might not be even necessary in certain cases, if the file is included implicitely by another header file
 * or by the test itself. It's safer though. 
 * 
 */
#include "ParallelTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "MonodomainPde.hpp"
#include "BidomainPde.hpp"


/**
 * Base class for cardiac problems; contains code generic to both mono- and bidomain.
 *
 * See tutorials for usage.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractCardiacProblem
{
    friend class TestBidomainWithBathAssembler;
    friend class TestCardiacSimulationArchiver;

    /** To save typing */
    typedef typename boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> >
        BccType;
    
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    
    /**
     * Save the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & mMeshFilename;
        archive & mpMesh;
        //archive & mAllocatedMemoryForMesh; // Mesh is deleted by AbstractCardiacPde
        archive & mNodesPerProcessorFilename;
        archive & mUseMatrixBasedRhsAssembly;
        archive & mWriteInfo;
        archive & mPrintOutput;
        archive & mNodesToOutput;
        //archive & mVoltageColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mExtraVariablesId; // Created by InitialiseWriter, called from Solve
        //archive & mTimeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mNodeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mpWriter; // Created by InitialiseWriter, called from Solve
        archive & mpCardiacPde;
        //archive & mpAssembler; // Only exists during calls to the Solve method
        bool has_solution = (mSolution != NULL);
        archive & has_solution;
        if (has_solution)
        {
            /// \todo: code for saving/loading mSolution is PROBLEM_DIM specific, move it into the save/load methods fo Mono and BidomainProblem                        
            std::string filename = ArchiveLocationInfo::GetArchiveDirectory() + "AbstractCardiacProblem_mSolution.vec";

            Hdf5DataWriter writer(*mpMesh->GetDistributedVectorFactory(), ArchiveLocationInfo::GetArchiveRelativePath(), "AbstractCardiacProblem_mSolution", false);
            writer.DefineFixedDimension(mpMesh->GetDistributedVectorFactory()->GetProblemSize());
            
            writer.DefineUnlimitedDimension("Time", "msec");
            int vm_col = writer.DefineVariable("Vm","mV");              
            
            if (PROBLEM_DIM==1)
            {
                writer.EndDefineMode();
                writer.PutUnlimitedVariable(0.0);           
                writer.PutVector(vm_col, mSolution);    
            }
            
            if (PROBLEM_DIM==2)
            {
                int phie_col = writer.DefineVariable("Phie","mV");
                writer.EndDefineMode();                           
                writer.PutUnlimitedVariable(0.0);           
                writer.PutStripedVector(vm_col, phie_col, mSolution);    
            }                
            
            writer.Close();            
            
        }
        archive & mCurrentTime;
        
        // Save boundary conditions
        SaveBoundaryConditions(archive, mpMesh, mpBoundaryConditionsContainer);
        SaveBoundaryConditions(archive, mpMesh, mpDefaultBoundaryConditionsContainer);
    }
    
    /**
     * Load the member variables.
     * 
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & mMeshFilename;
        archive & mpMesh;
        //archive & mAllocatedMemoryForMesh; // Will always be true after a load
        archive & mNodesPerProcessorFilename;
        archive & mUseMatrixBasedRhsAssembly;
        archive & mWriteInfo;
        archive & mPrintOutput;
        archive & mNodesToOutput;
        //archive & mVoltageColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mExtraVariablesId; // Created by InitialiseWriter, called from Solve
        //archive & mTimeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mNodeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mpWriter; // Created by InitialiseWriter, called from Solve
        archive & mpCardiacPde;
        //archive & mpAssembler; // Only exists during calls to the Solve method
        bool has_solution;
        archive & has_solution;
        if (has_solution)
        {
            /// \todo: code for saving/loading mSolution is PROBLEM_DIM specific, move it into the save/load methods fo Mono and BidomainProblem                        
            std::string filename = ArchiveLocationInfo::GetArchiveDirectory() + "AbstractCardiacProblem_mSolution.vec";

            mSolution = mpMesh->GetDistributedVectorFactory()->CreateVec(PROBLEM_DIM);
            DistributedVector mSolution_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(mSolution);                    

            Vec vm = mpMesh->GetDistributedVectorFactory()->CreateVec();
            Vec phie = mpMesh->GetDistributedVectorFactory()->CreateVec();

            Hdf5DataReader reader(ArchiveLocationInfo::GetArchiveRelativePath(), "AbstractCardiacProblem_mSolution", ArchiveLocationInfo::GetIsDirRelativeToChasteTestOutput());
            reader.GetVariableOverNodes(vm, "Vm", 0);

            if (PROBLEM_DIM==1)
            {
                //reader.Close(); // no need to call close explicitly, done in the destructor

                DistributedVector vm_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(vm);                    
                
                DistributedVector::Stripe mSolution_vm(mSolution_distri,0);
    
                for (DistributedVector::Iterator index = mSolution_distri.Begin();
                     index != mSolution_distri.End();
                     ++index)
                {
                    mSolution_vm[index] = vm_distri[index];
                }
            }
            
            if (PROBLEM_DIM==2)
            {
                reader.GetVariableOverNodes(phie, "Phie", 0);
                //reader.Close(); // no need to call close explicitly, done in the destructor
    
                DistributedVector vm_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(vm);                    
                DistributedVector phie_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(phie);
                
                DistributedVector::Stripe mSolution_vm(mSolution_distri,0);
                DistributedVector::Stripe mSolution_phie(mSolution_distri,1);
    
                for (DistributedVector::Iterator index = mSolution_distri.Begin();
                     index != mSolution_distri.End();
                     ++index)
                {
                    mSolution_vm[index] = vm_distri[index];
                    mSolution_phie[index] = phie_distri[index];
                }
            }

            mSolution_distri.Restore();
            
            VecDestroy(vm);
            VecDestroy(phie);            

        }
        archive & mCurrentTime;
        
        // Load boundary conditions
        mpBoundaryConditionsContainer = LoadBoundaryConditions(archive, mpMesh);
        mpDefaultBoundaryConditionsContainer = LoadBoundaryConditions(archive, mpMesh);
    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
    /**
     * Serialization helper method to save a boundary conditions container.
     *
     * @param archive  the archive to save to
     * @param pMesh  the mesh boundary conditions are defined on
     * @param pBcc  the container to save
     */
    template<class Archive>
    void SaveBoundaryConditions(Archive & archive,
                                AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                BccType pBcc) const
    {
        ///\todo #1169 need more tests of this!
        (*ProcessSpecificArchive<Archive>::Get()) & pBcc;
    }
    
    /**
     * Serialization helper method to load a boundary conditions container.
     *
     * @param archive  the archive to load from
     * @param pMesh  the mesh boundary conditions are to be defined on
     * @return  the loaded container
     */
    template<class Archive>
    BccType LoadBoundaryConditions(
            Archive & archive,
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
                                
    {
        // Load pointer from archive
        BccType p_bcc;
        (*ProcessSpecificArchive<Archive>::Get()) & p_bcc;
        
        // Fill in the conditions, if we have a container and it's not already full
        if (p_bcc)
        {
            p_bcc->LoadFromArchive(*ProcessSpecificArchive<Archive>::Get(), pMesh);
        }

        return p_bcc;
    }

protected:
    /** Meshes can be read from file or instantiated and passed directly to this
     *  class, this is for the former */
    std::string mMeshFilename;

    /** If this is set, the nodes for each processor are read */
    std::string mNodesPerProcessorFilename;

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

    /** If only outputing voltage for selected nodes, which nodes to output at */
    std::vector<unsigned> mNodesToOutput;

    /** Used by the writer */
    unsigned mVoltageColumnId;
    /** List of extra variables to be written to HDF5 file */
    std::vector<unsigned> mExtraVariablesId;
    /** Used by the writer */
    unsigned mTimeColumnId;
    /** Used by the writer */
    unsigned mNodeColumnId;

    /** The monodomain or bidomain pde */
    AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>* mpCardiacPde;

    /** Boundary conditions container used in the simulation */
    BccType mpBoundaryConditionsContainer;
    /** It is convenient to also have a separate variable for default (zero-Neumann) boundary conditions */
    BccType mpDefaultBoundaryConditionsContainer;
    /** The PDE solver */
    AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpAssembler;
    /** The cell factory creates the cells for each node */
    AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* mpCellFactory;
    /** The mesh. Can either by passed in, or the mesh filename can be set */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /** The current solution vector, of the form [V_0 .. V_N ] for monodomain and
     *  [V_0 phi_0 .. V_N phi_N] for bidomain */
    Vec mSolution;
    
    /**
     * The current simulation time.
     * 
     * This is used to be able to restart simulations at a point other than time zero,
     * either because of repeated calls to Solve (with increased simulation duration)
     * or because of restarting from a checkpoint.
     */
    double mCurrentTime;

    /**
     * Subclasses must override this method to create a PDE object of the appropriate type.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>* CreateCardiacPde() =0;

    /**
     * Subclasses must override this method to create a suitable assembler object.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* CreateAssembler() =0;

public:
    /**
     * The object to use to write results to disk.
     * 
     * This (and things in MonodomainProblem) being public are hacks for
     * CardiacElectroMechanicsProblem to work.
     * 
     * \todo CardiacElectroMechanicsProblem should be a friend, but not sure
     * how to get friends to work when both friends are templated and abstract.
     */
    Hdf5DataWriter* mpWriter;

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    AbstractCardiacProblem(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory);
    
    /**
     * Constructor used by archiving.
     */
    AbstractCardiacProblem();
    
    /**
     *  Destructor
     */
    virtual ~AbstractCardiacProblem();

    /**
     * Initialise the system, once parameters have been set up.
     * 
     * Must be called before first calling Solve().  If loading from a checkpoint,
     * do NOT call this method, as it can also be used to reset the problem to
     * perform another simulation from time 0.
     */
    void Initialise();

    /**
     *  Set a file from which the nodes for each processor are read
     * 
     * @param rFilename
     */
    void SetNodesPerProcessorFilename(const std::string& rFilename);

    /**
     *  Set the boundary conditions container.
     *  @param pBcc is a pointer to a boundary conditions container
     */
    void SetBoundaryConditionsContainer(BccType pBcc);

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
     * This only needs to be called if a mesh filename has not been set.
     * 
     * @param pMesh  the mesh object to use
     */
    void SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     *  Set whether the simulation will generate results files.
     * 
     * @param rPrintOutput
     */
    void PrintOutput(bool rPrintOutput);

    /**
     *  Set whether extra info will be written to stdout during computation.
     * 
     * @param writeInfo
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
    
    /**
     * Get the solution vector, wrapped in a DistributedVector.
     * 
     * See also GetSolution.
     */
    DistributedVector GetSolutionDistributedVector();

    /**
     * @return the mesh used
     */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM> & rGetMesh();

    /**
     * @return the cardiac PDE used
     */
    AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>* GetPde();

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

    /**
     * Write informative details about the progress of the simulation to standard output.
     * 
     * Implemented only in subclasses.
     * 
     * @param time  the current time
     */
    virtual void WriteInfo(double time)=0;

    /**
     * Define what variables are written to the primary results file.
     * @param extending  whether we are extending an existing results file
     */
    virtual void DefineWriterColumns(bool extending);
    
    /**
     * Define the user specified variables to be written to the primary results file
     * @param extending  whether we are extending an existing results file
     */
    void DefineExtraVariablesWriterColumns(bool extending);

    /**
     * Write one timestep of output data to the primary results file.
     * 
     * @param time  the current time
     * @param voltageVec  the solution vector to write
     */
    virtual void WriteOneStep(double time, Vec voltageVec) = 0;
    
    /**
     * Write one timestep of output data for the extra variables to the primary results file.
     */
    void WriteExtraVariablesOneStep();

    /**
     * It creates and initialises the hdf writer from the Hdf5DataWriter class.
     * It passes the output directory and file name to it.
     * It is called by Solve(), if the output needs to be generated.
     */
    void InitialiseWriter();

    /**
     * Specifies which nodes in the mesh to output.
     * 
     * @param rNodesToOutput is a reference to a vector with the indexes of the nodes
     * where the output is desired.
     * If empty, the output will be for all the nodes in the mesh.
     */
    void SetOutputNodes(std::vector<unsigned> & rNodesToOutput);

    /**
     * Create and return a data reader configured to read the results we've been outputting.
     */
    Hdf5DataReader GetDataReader();

    /**
     *  Whether to use matrix-based RHS assembly or not.
     * 
     * @param usematrix
     */
    void UseMatrixBasedRhsAssembly(bool usematrix=true);

    /**
     *  Called at end of each time step in the main time-loop in
     *  Solve(). Empty implementation but can be overloaded by child
     *  classes.
     * 
     * @param time  the current time
     */
    virtual void OnEndOfTimestep(double time)
    {}
    
    /**
     *  Tells the problem class to archive the linear system after every call to Solve()
     * 
     * @param  archive set true to archive the LinearSystem object at the end
     */
    void SetArchiveLinearSystemObject(bool archive=true);
    
    /**
     * Used when loading a set of archives written by a parallel simulation onto a single process.
     * Loads data from the given process-specific archive (written by a non-master process) and
     * merges it into our data.
     * 
     * @param archive  the archive to load.
     * 
     * \note The process-specific archives currently contain the following data.  If the layout changes,
     * then this method will need to be altered, since it hard-codes knowledge of the order in
     * which things are archived.
     * 
     *  -# (via #mpMesh) DistributedVectorFactory*
     *  -# (via #mpCardiacPde LoadCardiacCells) DistributedVectorFactory*
     *  -# (via #mpCardiacPde LoadCardiacCells) number_of_cells and sequence of AbstractCardiacCell*
     *  -# (via #mpCardiacPde) DistributedVectorFactory*
     *  -# #mpBoundaryConditionsContainer
     *  -# #mpDefaultBoundaryConditionsContainer
     */
    template<class Archive>
    void LoadExtraArchive(Archive & archive);
};

TEMPLATED_CLASS_IS_ABSTRACT_3_UNSIGNED(AbstractCardiacProblem);

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
template<class Archive>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LoadExtraArchive(Archive & archive)
{
    // The vector factory must be loaded, but isn't needed for anything.
    DistributedVectorFactory* p_mesh_factory;
    archive >> p_mesh_factory;
    
    {
        // The cells collection vector factory
        DistributedVectorFactory* p_cells_factory;
        archive >> p_cells_factory;
        assert(p_cells_factory == p_mesh_factory); // Paranoia...
    }

    // The cardiac cells
    std::vector<AbstractCardiacCell*> cells;
    unsigned num_cells;
    archive >> num_cells;
    cells.reserve(num_cells);
    for (unsigned i=0; i<num_cells; i++)
    {
        AbstractCardiacCell* p_cell;
        archive >> p_cell;
        cells.push_back(p_cell);
    }
    mpCardiacPde->ExtendCells(cells);

    {
        DistributedVectorFactory* p_pde_factory;
        archive >> p_pde_factory;
        assert(p_pde_factory == p_mesh_factory); // Paranoia...
    }
    // We no longer need this vector factory, since we already have our own.
    delete p_mesh_factory;

    // The boundary conditions
    BccType p_bcc;
    archive >> p_bcc;
    if (p_bcc)
    {
        if (!mpBoundaryConditionsContainer)
        {
            mpBoundaryConditionsContainer = p_bcc;
        }
        mpBoundaryConditionsContainer->MergeFromArchive(archive, mpMesh);
    }
    BccType p_default_bcc;
    archive >> p_default_bcc;
    if (p_default_bcc)
    {
        // This always holds, so we never need the commented code below.
        assert(p_bcc == p_default_bcc);
//        if (!mpDefaultBoundaryConditionsContainer)
//        {
//            mpDefaultBoundaryConditionsContainer = p_bcc;
//        }
//        mpDefaultBoundaryConditionsContainer->MergeFromArchive(archive, mpMesh);
    }
}


#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
