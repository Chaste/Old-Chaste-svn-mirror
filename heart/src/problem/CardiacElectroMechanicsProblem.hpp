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


#ifndef CARDIACELECTROMECHANICSPROBLEM_HPP_
#define CARDIACELECTROMECHANICSPROBLEM_HPP_

#include <vector>
#include <string>
#include "UblasIncludes.hpp"

#include "AbstractCardiacCellFactory.hpp"
#include "MonodomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "AbstractOdeBasedContractionModel.hpp"
#include "AbstractCardiacMechanicsSolver.hpp"
#include "AbstractCardiacMechanicsSolverInterface.hpp"
#include "FineCoarseMeshPair.hpp"
#include "AbstractConductivityModifier.hpp"

/**
 *  CardiacElectroMechanicsProblem
 *
 *  For solving full electro-mechanical problems.
 *
 *  Solves a monodomain problem (diffusion plus cell models) on a (fine) electrics
 *  mesh, and a mechanics problem (finite elasticity plus contraction model) on a coarse
 *  mesh. An implicit scheme (Jon Whiteley's algorithm) be be used.
 *
 *  For solving problems on regular grids use CardiacElectroMechProbRegularGeom
 *
 *  The implicit algorithm:
 *
 *  Store the position in the electrics mesh of each quad point in the mechanics mesh
 *  For every time:
 *    Solve the monodomain problem (ie integrate ODEs, solve PDE)
 *    Get intracellular [Ca] at each electrics node and interpolate on each mechanics quad point
 *    Set [Ca] on each contraction model (one for each point)
 *    Solve static finite elasticity problem implicity
 *       - guess solution
 *       - this gives the fibre stretch and stretch rate to be set on contraction models
 *       - integrate contraction models (implicitly if NHS) to get for active tension
 *       - use this active tension in computing the stress for that guess of the deformation
 *  end
 */
template<unsigned DIM>
class CardiacElectroMechanicsProblem
    : public AbstractConductivityModifier<DIM,DIM> // this only inherits from this class so it can be passed to the tissue to
                                                   // allow deformation-based altering of the conductivity
{

friend class TestCardiacElectroMechanicsProblem;

protected :
    /** Contraction model (from enumeration) */
    ContractionModel mContractionModel;
    /** The cardiac problem class */
    MonodomainProblem<DIM>* mpMonodomainProblem;

///\todo #1756
    /** The mechanics solver - a pointer to the part that sees the cardiac mechanics interface bit */
    AbstractCardiacMechanicsSolverInterface<DIM>* mpCardiacMechSolver;
    /** The mechanics solver - a pointer to the part that sees the solid mechanics solver */
    AbstractNonlinearElasticitySolver<DIM>* mpMechanicsSolver;

    /** End time. The start time is assumed to be 0.0 */
    double mEndTime;
    /** The electrics timestep. */
    double mElectricsTimeStep;
    /** The mechanics timestep. Needs to be a multiple of the electrics timestep */
    double mMechanicsTimeStep;
    /** The number of electrics timesteps per mechanics timestep */
    unsigned mNumElecTimestepsPerMechTimestep;
    /** Timestep to use when solving contraction models */
    double mContractionModelOdeTimeStep;

    /** The mesh for the electrics */
    TetrahedralMesh<DIM,DIM>* mpElectricsMesh;
    /** The mesh for the mechanics */
    QuadraticMesh<DIM>* mpMechanicsMesh;

    /**
     *  Object containing body force and boundary conditions info for
     *  solid mechanics problems.
     */
    SolidMechanicsProblemDefinition<DIM>* mpProblemDefinition;

    /** Class wrapping both meshes, useful for transferring information */
    FineCoarseMeshPair<DIM>* mpMeshPair;

    /** Output directory, relative to TEST_OUTPUT */
    std::string mOutputDirectory;
    /** Deformation output-sub-directory */
    std::string mDeformationOutputDirectory;
    /** Whether to write any output */
    bool mWriteOutput;
    /** Whether to not write out voltages */
    bool mNoElectricsOutput;

    /** when to write output */
    static const int WRITE_EVERY_NTH_TIME = 1; //hardcoded for the time being ///\todo, allow user to set this

    /** Whether any location has been set to be watched (lots of output for that location */
    bool mIsWatchedLocation;
    /** The watched location if there is one */
    c_vector<double,DIM> mWatchedLocation;
    /** The node in the electrics mesh corresponding to the watched location */
    unsigned mWatchedElectricsNodeIndex;
    /** The node in the mechanics mesh corresponding to the watched location */
    unsigned mWatchedMechanicsNodeIndex;
    /** File where watched location info is written */
    out_stream mpWatchedLocationFile;

    /** Nodes for which the deformation is fixed to zero */
    std::vector<unsigned> mFixedNodes;
    /** .ortho/.orthoquad file from which to read element-wise/quadpoint-wise fibre-sheet-normal-directions */
    std::string mFibreSheetDirectionsFile;
    /** Whether the mFibreSheetDirectionsFile file gives the fibre info at each element, or each quadrature point */
    bool mFibreSheetDirectionsDefinedPerQuadraturePoint;

    /** A vector of stretches (in the fibre direction), one for each element in the mechanics mesh */
    std::vector<double> mStretchesForEachMechanicsElement;
    /** A vector of deformation gradients (each entry a matrix), one for each element in the mechanics mesh */
    std::vector<c_matrix<double,DIM,DIM> > mDeformationGradientsForEachMechanicsElement;

    /** A pair of (element_index, deformed_conductivity) for the last element on which the deformed conductivity
     *  sigma_def = F^{-1} sigma_undef F^{-T} has been computed. Used in rGetModifiedConductivityTensor().
     */
    std::pair<unsigned, c_matrix<double,DIM,DIM> > mLastModifiedConductivity;

    /** Whether the deformation is always to alter the conductivities (ie one part of MEF) */
    bool mConductivityAffectedByDeformationMef;

    /** Whether the deformation is always to affect the cell models (for example, for use in stretch-activated channels)
     *  (ie one part of MEF)
     */
    bool mCellModelsAffectedByDeformationMef;

    /** The material law to be used. Defaults to NashHunterPoleZero if an incompressible problem is being solved */
    AbstractMaterialLaw<DIM>* mpMaterialLaw;

    /** Whether a material law was passed in or the default used */
    bool mAllocatedMaterialLawMemory;


    /**
     *  Determine which node is closest to the watched location
     */
    void DetermineWatchedNodes();


    /**
     *  Write info (x, y, [z], V) for the watched node.
     *
     * @param time  Time-step now, to write out
     * @param voltage  Vm vector (this is Monodomain)
     */
    void WriteWatchedLocationData(double time, Vec voltage);


//// #1245
//    std::vector<BoundaryElement<DIM-1,DIM>*>* mpImpactRegion;
//    std::vector<c_vector<double,DIM> > mImpactTractions;
//    void ApplyImpactTractions(double time);

public :

    /**
     * Constructor.
     * @param contractionModel contraction model (see the enum "ContractionModel" for the options).
     * @param pElectricsMesh  Mesh on which to solve electrics (Monodomain)
     * @param pMechanicsMesh  Mesh (2nd order) on which to solve mechanics
     * @param fixedMechanicsNodes  Indices of those nodes which a pinned in space
     * @param pCellFactory factory to use to create cells
     * @param endTime the end time to use
     * @param electricsPdeTimeStep timestep used in solving for the electrical activity
     * @param mechanicsSolveTimestep how often the mechanics is solved for (should be a multiple of electrics PDE timestep)
     * @param contractionModelOdeTimeStep Step size for contraction model (of active tension in cardiac cells) being used.
     * @param outputDirectory the output directory
     */
    CardiacElectroMechanicsProblem(ContractionModel contractionModel,
                                   TetrahedralMesh<DIM,DIM>* pElectricsMesh,
                                   QuadraticMesh<DIM>* pMechanicsMesh,
                                   std::vector<unsigned> fixedMechanicsNodes,
                                   AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   double endTime,
                                   double electricsPdeTimeStep,
                                   double mechanicsSolveTimestep,
                                   double contractionModelOdeTimeStep,
                                   std::string outputDirectory);

    /**
     *  Delete allocated memory and close the watched location file
     *
     *  NOTE if SetWatchedLocation but not Initialise has been
     *  called, mpWatchedLocationFile will be uninitialised and
     *  using it will cause a seg fault. Hence the mpMechanicsMesh!=NULL
     *  it is true if Initialise has been called.
     */
    virtual ~CardiacElectroMechanicsProblem();

    /**
     *  Initialise the class. Initialises the MonodomainProblem
     *  and sets up the electrics mesh to mechanics mesh data.
     */
    void Initialise();

    /**
     *  Solve the electromechanics problem
     */
    void Solve();

    /**
     *  Short helper function - the max of a std::vector
     *  @param vec a vector of doubles
     */
    double Max(std::vector<double>& vec);

    /** Call to not write out voltages */
    void SetNoElectricsOutput();

    /**
     *  Set a location to be watched - for which lots of output
     *  is given. Should correspond to nodes in both meshes.
     *
     *  The watched file will have rows that look like:
     *  time x_pos y_pos [z_pos] voltage Ca_i_conc.
     *
     *  NOTE: for the Calcium - assumes LUO_RUDY IS USED
     *  @param watchedLocation  location (x,y,z) in space.  Watched node is the closest to this point.
     */
    void SetWatchedPosition(c_vector<double,DIM> watchedLocation);

    /**
     *  Set a variable fibre-sheet-normal direction (matrices), from file.
     *  If the second parameter is false, there should be one fibre-sheet definition for each element; otherwise
     *  there should be one fibre-sheet definition for each *quadrature point* in the mesh.
     *  In the first case, the file should be a .ortho file (ie each line has the fibre dir, sheet dir, normal dir 
     *  for that element), in the second it should have .orthoquad as the format.
     * 
     *  @param orthoFile the file containing the fibre/sheet directions
     *  @param definedPerQuadPoint whether the fibre-sheet definitions are for each quadrature point in the mesh
     *   (if not, one for each element is assumed).
     */
    void SetVariableFibreSheetDirectionsFile(std::string orthoFile, bool definedPerQuadPoint);

    /** @return the current deformed position of the nodes */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition();

    /**
     *  By default (at the moment), the deformation does not affect the electrophysiology in any way. 
     *  Call this to allow it to, then
     *   (i)  The stretch will be passed back to the cell models for use stretch-activated channels etc
     *   (ii) The deformation to alter the conductivity
     * 
//EMTODO - check these
     *  Two things to note:
     *   (i) this can't be called if fibre-sheet directions have been defined from file for each quadrature
     *  point (as opposed to each mechanics element) - this is because if the stretch is to be passed back to
     *  the electric mesh nodes, the fibre direction has to be defined at those nodes
     *   (ii) currently the set-up stage (computing mechanics mesh elements and weights for electrics mesh 
     *  nodes) is inefficiently implemented - setup will be very slow for big meshes
     */
    void UseMechanoElectricFeedback()
    {
    	mConductivityAffectedByDeformationMef = true;
    	mCellModelsAffectedByDeformationMef = true;
    }

    /**
     *  The implementation of the pure method defined in the base class AbstractConductivityModifier. The tissue class will
     *  call this method.
     *  @param elementIndex Index of current element
     *  @param rOriginalConductivity Reference to the original (for example, undeformed) conductivity tensor
     *  @return Reference to a modified conductivity tensor.
     */
    c_matrix<double,DIM,DIM>& rGetModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,DIM,DIM>& rOriginalConductivity);

//// #1245
//    void SetImpactRegion(std::vector<BoundaryElement<DIM-1,DIM>*>& rImpactRegion);

    /**
     * Called in Solve() before the time loop
     */
    virtual void PrepareForSolve(){}

    /**
     * Called in Solve() at the end of every time step
     * @param counter time step
     */
    virtual void OnEndOfTimeStep(unsigned counter){}

    /**
     *  Get a pointer to the object which contains the mechanics problem definition (which includes
     *  boundary conditions, both fixed node or traction). The user can call this to say
     *  then set traction boundary conditions for the mechanics problem by calling methods on
     *  this object.
     *
     *  Needs to be called after Initialise().
     */
    SolidMechanicsProblemDefinition<DIM>* GetSolidMechanicsProblemDefinition()
    {
        assert(mpProblemDefinition);
        return mpProblemDefinition;
    }

    /** Set the material law. If this isn't called the default material law will be used.
     *  Only call before Initialise() and before Solve()
     *  @param pMaterialLaw the material law
     */
    void SetMaterialLaw(AbstractMaterialLaw<DIM>* pMaterialLaw)
    {
        assert(pMaterialLaw);
        assert(mpMaterialLaw == NULL); // only call SetMaterialLaw() before both Initialise() or Solve()
        mpMaterialLaw = pMaterialLaw;
        mAllocatedMaterialLawMemory = false;
    }
};



#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
