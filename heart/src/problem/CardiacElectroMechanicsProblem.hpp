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

//// TODO: 
//     more tests 
//     use direct solver
//     go through and tidy/refactor, perhaps make elements and weights safer


#ifndef CARDIACELECTROMECHANICSPROBLEM_HPP_
#define CARDIACELECTROMECHANICSPROBLEM_HPP_

#include <vector>
#include <string>
#include "UblasIncludes.hpp"

#include "AbstractCardiacCellFactory.hpp"
#include "MonodomainProblem.hpp"
#include "ImplicitCardiacMechanicsAssembler.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"

// if including Cinv in monobidomain equations
//#include "NodewiseData.hpp"


/* todos:
 *
 * think about architecture (of AbstractCardiacProblem) when this is done properly..
 */



/**
 *  At the beginning of a two mesh simulation we need to figure out and store
 *  which (electrics-mesh) element each (mechanics-mesh) gauss point is in, and
 *  what the weight of that gauss point for that particular element is. This struct
 *  just contains this two pieces of data
 */
template<unsigned DIM>
struct ElementAndWeights
{
    unsigned ElementNum;
    c_vector<double,DIM+1> Weights;
};



/**
 *  CardiacElectroMechanicsProblem
 *
 *  For solving full electro-mechanical problems. 
 *
 *  Solves a monodomain problem (diffusion plus cell models) on a (fine) electrics 
 *  mesh, and a mechanics problem (finite elasticity plus NHS cell models) on a coarse 
 *  mesh. An implicit scheme (Jon Whiteley's algorithm) be be used.
 *
 *  At the moment just solves on the unit square. The spatial stepsize for the
 *  electrics is fixed (96 by 96), and the displacement boundary conditions are
 *  zero displacement on X=0.
 *
 *  The implicit algorithm:
 *
 *  Store the position in the electrics mesh of each quad point in the mechanics mesh
 *  For every time:
 *    Solve the monodomain problem (ie integrate ODEs, solve PDE)
 *    Get intracellular [Ca] at each electrics node and interpolate on each mechanics quad point
 *    Set [Ca] on each NHS model (one for each point)
 *    Solve static finite elasticity problem implicity
 *       - guess solution
 *       - this gives the fibre stretch and stretch rate to be set on NHS models
 *       - integrate NHS models implicity for active tension
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
    /*< The cardiac problem class */
    MonodomainProblem<DIM>* mpMonodomainProblem;
    /*< The mechanics assembler */
    ImplicitCardiacMechanicsAssembler<DIM>* mpCardiacMechAssembler;

    /*< End time. The start time is assumed to be 0.0 */
    double mEndTime;
    /*< The electrics timestep. */
    double mElectricsTimeStep;
    /*< The mechanics timestep. Needs to be a multiple of the electrics timestep */
    double mMechanicsTimeStep;
    /*< The number of electrics timesteps per mechanics timestep */
    unsigned mNumElecTimestepsPerMechTimestep;
    /*< Timestep to use when solving NHS models (for implicit version)*/
    double mNhsOdeTimeStep;

    /*< The mesh for the electrics */
    TetrahedralMesh<DIM,DIM>* mpElectricsMesh;
    /*< The mesh for the mechanics */
    QuadraticMesh<DIM>* mpMechanicsMesh;

    /**
     *  The (electrics-mesh) element numbers saying which element each
     *  (mechanics-mesh) gauss point is in, and the weight of that gauss point
     *  for that particular element.
     */
    std::vector<ElementAndWeights<DIM> > mElementAndWeightsForQuadPoints;

    /*< Output directory, relative to TEST_OUTPUT */
    std::string mOutputDirectory;
    /*< Deformation output-sub-directory */
    std::string mDeformationOutputDirectory;
    /*< Whether to write any output */
    bool mWriteOutput;
    /** Whether to not write out voltages */
    bool mNoElectricsOutput;
    /*< when to write output */
    const static int WRITE_EVERY_NTH_TIME = 1;

    /**
     *  Whether to use a direct solver when solving linear system. Should
     *  definitely be used if UMFPACK is installed.
     */
    bool mUseDirectLinearSolver;

    /*< Whether any location has been set to be watched (lots of output for that location */
    bool mIsWatchedLocation;
    /*< The watched location if there is one */
    c_vector<double,DIM> mWatchedLocation;
    /*< The node in the electrics mesh corresponding to the watched location */
    unsigned mWatchedElectricsNodeIndex;
    /*< The node in the mechanics mesh corresponding to the watched location */
    unsigned mWatchedMechanicsNodeIndex;
    /*< File where watched location info is written */
    out_stream mpWatchedLocationFile;

    /*< Nodes for which the deformation is fixed to zero */
    std::vector<unsigned> mFixedNodes;

    /** 
     *  Determine which node is closest to the watched location
     */
    void DetermineWatchedNodes();


    /**
     *  Write info (x, y, V, and Ca) for the watched node. Note: the Ca is written, 
     *  but this ASSUMES LUO-RUDY IS USED
     */
    void WriteWatchedLocationData(double time, Vec voltage);

public :
    /**
     *  Constructor
     */
    CardiacElectroMechanicsProblem(TetrahedralMesh<DIM,DIM>* pElectricsMesh,
                                   QuadraticMesh<DIM>* pMechanicsMesh,
                                   std::vector<unsigned> fixedMechanicsNodes,
                                   AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   double endTime,
                                   unsigned numElecTimeStepsPerMechTimestep,
                                   double nhsOdeTimeStep,
                                   std::string outputDirectory);

    /**
     *  Delete allocated memory and close the watched location file
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

    // short helper function - the max of a std::vec. this is in the wrong place
    double Max(std::vector<double>& vec);

    /** Call to not write out voltages */
    void SetNoElectricsOutput();

//    /** Use the direct solver when solving linear systems in the
//     *  mechanics. DEFINITELY should be used in experimental work.
//     */
//    void UseDirectLinearSolver();

    /**
     *  Set a location to be watched - for which lots of output
     *  is given. Should correspond to nodes in both meshes.
     *
     *  The watched file will have rows that look like:
     *  time x_pos y_pos [z_pos] voltage Ca_i_conc.
     *  
     *  NOTE: for the Calcium - assumes LUO_RUDY IS USED
     */
    void SetWatchedPosition(c_vector<double,DIM> watchedLocation);
    
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition();
};



#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
