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
#include "AbstractCardiacMechanicsAssembler.hpp"

// if including Cinv in monobidomain equations
//#include "NodewiseData.hpp"


// EMTODO: make elements and weights safer?



/**
 *  At the beginning of a two mesh simulation we need to figure out and store
 *  which (electrics-mesh) element each (mechanics-mesh) gauss point is in, and
 *  what the weight of that gauss point for that particular element is. This struct
 *  just contains this two pieces of data
 */
template<unsigned DIM>
struct ElementAndWeights
{
    unsigned ElementNum; /**< Which element*/
    c_vector<double,DIM+1> Weights; /**<Gauss weights for this element*/
};



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
 *
 *  Note: invC is not used in the monodomain equations (code added but commented
 *  out) we have shown that this does not affect the mechanics results (might
 *  affect the electrics).
 */
template<unsigned DIM>
class CardiacElectroMechanicsProblem
{

friend class TestCardiacElectroMechanicsProblem;

protected :
    /** Contraction model (from enumeration) */
    ContractionModel mContractionModel;
    /** The cardiac problem class */
    MonodomainProblem<DIM>* mpMonodomainProblem;
    /** The mechanics assembler */
    AbstractCardiacMechanicsAssembler<DIM>* mpCardiacMechAssembler;

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
     *  The (electrics-mesh) element numbers saying which element each
     *  (mechanics-mesh) gauss point is in, and the weight of that gauss point
     *  for that particular element.
     */
    std::vector<ElementAndWeights<DIM> > mElementAndWeightsForQuadPoints;

    /** Output directory, relative to TEST_OUTPUT */
    std::string mOutputDirectory;
    /** Deformation output-sub-directory */
    std::string mDeformationOutputDirectory;
    /** Whether to write any output */
    bool mWriteOutput;
    /** Whether to not write out voltages */
    bool mNoElectricsOutput;

    /** when to write output */
    const static int WRITE_EVERY_NTH_TIME = 1; //hardcoded for the time being ///\todo, allow user to set this

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
	/** .ortho file from which to read element-wise fibre-sheet-normal-directions */
    std::string mFibreSheetDirectionsFile;

    /**
     *  Determine which node is closest to the watched location
     */
    void DetermineWatchedNodes();


//EMTODO
    /**
     *  Write info (x, y, V, and Ca) for the watched node. Note: the Ca is written,
     *  but this ASSUMES LUO-RUDY IS USED
     * 
     * @param time  Time-step now, to write out
     * @param voltage  Vm vector (this is Monodomain)
     */
    void WriteWatchedLocationData(double time, Vec voltage);

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
     * @param numElecTimeStepsPerMechTimestep number of electrics timesteps to be used in each mechanics solve
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
                                   unsigned numElecTimeStepsPerMechTimestep,
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
     *  Initialise the class. Calls ConstructMeshes() and
     *  ConstructMechanicsAssembler(). Initialises the MonodomainProblem
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
     *  Set fibre/sheet directions for each element from a file. 
     *  The file should be a .ortho file (ie each line has the fibre dir, sheet dir, normal dir for that element).
     *  The number of elements must match the number in the MECHANICS mesh!
     *  @param orthoFile the file containing the fibre/sheet directions
     */
    void SetVariableFibreSheetDirectionsFile(std::string orthoFile);

    /** @return the current deformed position of the nodes */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition();
};



#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
