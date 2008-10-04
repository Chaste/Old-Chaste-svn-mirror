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
//     more tests (of this and implicit assembler)
//     run for longer time
//     compare times
//     use direct solver
//     go through and tidy/refactor, perhaps make elements and weights safer


#ifndef CARDIACELECTROMECHANICSPROBLEM2_HPP_
#define CARDIACELECTROMECHANICSPROBLEM2_HPP_

#include "MonodomainProblem.hpp"
#include "ImplicitCardiacMechanicsAssembler2.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "NhsCellularMechanicsOdeSystem.hpp"
#include "TrianglesMeshWriter.hpp"
#include "NonlinearElasticityTools.hpp"
#include "QuadraturePointsGroup.hpp"

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
class CardiacElectroMechanicsProblem2
{

friend class TestCardiacElectroMechanicsProblem2;

protected :
    /*< The cardiac problem class */
    MonodomainProblem<DIM>* mpMonodomainProblem;
    /*< The mechanics assembler */
    ImplicitCardiacMechanicsAssembler2<DIM>* mpCardiacMechAssembler;

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
    ConformingTetrahedralMesh<DIM,DIM>* mpElectricsMesh;
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
    /*< Number of nodes per dimension in the mechanics mesh (taken in in constructor) */
    unsigned mNumElementsPerDimInMechanicsMesh;

    /**
     *  Construct the two meshes
     */
    void ConstructMeshes()
    {
        double width = 1.0;

        // create electrics mesh
        mpElectricsMesh = new ConformingTetrahedralMesh<DIM,DIM>();

        unsigned num_elem = 96;
        mpElectricsMesh->ConstructRectangularMesh(num_elem,num_elem);
        mpElectricsMesh->Scale(width/num_elem,width/num_elem);

        // create mechanics mesh
        assert(DIM==2); // the below assumes DIM==2 
        mpMechanicsMesh = new QuadraticMesh<DIM>(1.0,1.0,mNumElementsPerDimInMechanicsMesh,mNumElementsPerDimInMechanicsMesh);
        LOG(2, "Width of meshes is " << width);
        LOG(2, "Num nodes in electrical and mechanical meshes are: " << mpElectricsMesh->GetNumNodes() << ", " << mpMechanicsMesh->GetNumNodes() << "\n");
    }        
        
    
    /**
     *  Construct the mechanics assembler (the left-hand edge being fixed
     *  is hardcoded here, and the default material law is used
     */
    void ConstructMechanicsAssembler()
    {
        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<mpMechanicsMesh->GetNumNodes(); i++)
        {
            if( fabs(mpMechanicsMesh->GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }
        mpCardiacMechAssembler = new ImplicitCardiacMechanicsAssembler2<DIM>(mpMechanicsMesh,mDeformationOutputDirectory,fixed_nodes);
    }

    /** 
     *  Determine which node is closest to the watched location
     */
    void DetermineWatchedNodes()
    {
        assert(mIsWatchedLocation);

        // find the nearest electrics mesh node
        double min_dist = DBL_MAX;
        unsigned node_index = UNSIGNED_UNSET;
        for(unsigned i=0; i<mpElectricsMesh->GetNumNodes(); i++)
        {
            double dist = norm_2(mWatchedLocation - mpElectricsMesh->GetNode(i)->rGetLocation());
            if(dist < min_dist)
            {
                min_dist = dist;
                node_index = i;
            }
        }

        // set up watched node, if close enough
        assert(node_index != UNSIGNED_UNSET); // should def have found something
        c_vector<double,DIM> pos = mpElectricsMesh->GetNode(node_index)->rGetLocation();

        if(min_dist > 1e-8)
        {
            std::cout << "ERROR: Could not find an electrics node very close to requested watched location - "
                      << "min distance was " << min_dist << " for node " << node_index
                      << " at location " << pos << std::flush;;

            //// the following causes a seg fault for some reason (!!???!!!)
            // EXCEPTION("Could not find an electrics node very close to requested watched location");
            assert(0);
        }
        else
        {
            LOG_AND_COUT(2,"Chosen electrics node "<<node_index<<" at location " << pos << " to be watched");
            mWatchedElectricsNodeIndex = node_index;
        }

        // find nearest mechanics mesh
        min_dist = DBL_MAX;
        node_index = UNSIGNED_UNSET;
        c_vector<double,DIM> pos_at_min;

        for(unsigned i=0; i<mpMechanicsMesh->GetNumNodes(); i++)
        {
            c_vector<double,DIM> position = mpMechanicsMesh->GetNode(i)->rGetLocation();

            double dist = norm_2(position-mWatchedLocation);

            if(dist < min_dist)
            {
                min_dist = dist;
                node_index = i;
                pos_at_min = position;
            }
        }

        // set up watched node, if close enough
        assert(node_index != UNSIGNED_UNSET); // should def have found something

        if(min_dist > 1e-8)
        {
            std::cout << "ERROR: Could not find a mechanics node very close to requested watched location - "
                      << "min distance was " << min_dist << " for node " << node_index
                      << " at location " << pos_at_min;

            //// the following causes a seg fault for some reason (!!???!!!)
            //EXCEPTION("Could not find a mechanics node very close to requested watched location");
            assert(0);
        }
        else
        {
            LOG_AND_COUT(2,"Chosen electrics node "<<node_index<<" at location " << pos << " to be watched");
            mWatchedMechanicsNodeIndex = node_index;
        }

        OutputFileHandler handler(mOutputDirectory + "/watched/");
        mpWatchedLocationFile = handler.OpenOutputFile("data_NBassumesLr91.txt");
    }


    /**
     *  Write info (x, y, V, and Ca) for the watched node. Note: the Ca is written, 
     *  but this ASSUMES LUO-RUDY IS USED
     */
    void WriteWatchedLocationData(double time, Vec voltage)
    {
        assert(mIsWatchedLocation);

        std::vector<c_vector<double,DIM> >& deformed_position = mpCardiacMechAssembler->rGetDeformedPosition();

        ///\todo Rather inefficient
        ReplicatableVector voltage_replicated(voltage);
        double V=voltage_replicated[mWatchedElectricsNodeIndex];

///\todo: 
// NOTE!!! HARDCODED state variable index - assumes Lr91. Hierarchy not set up yet.
        double Ca = mpMonodomainProblem->GetMonodomainPde()->GetCardiacCell(mWatchedElectricsNodeIndex)->rGetStateVariables()[3];

        *mpWatchedLocationFile << time << " ";
        for(unsigned i=0; i<DIM; i++)
        {
            *mpWatchedLocationFile << deformed_position[mWatchedMechanicsNodeIndex](i) << " ";
        }
        *mpWatchedLocationFile << V <<  " " << Ca << "\n";
        mpWatchedLocationFile->flush();
    }

public :
    /**
     *  Constructor
     */
    CardiacElectroMechanicsProblem2(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                    double endTime,
                                    unsigned numElementsPerDimInMechanicsMesh,
                                    unsigned numElecTimeStepsPerMechTimestep,
                                    double nhsOdeTimeStep,
                                    std::string outputDirectory = "")
    {
        // create the monodomain problem. Note the we use this to set up the cells,
        // get an initial condition (voltage) vector, and get an assembler. We won't
        // ever call solve on the MonodomainProblem
        assert(pCellFactory != NULL);
        mpMonodomainProblem = new MonodomainProblem<DIM>(pCellFactory);

        // save time infomation
        assert(endTime > 0);
        mEndTime = endTime;
        mElectricsTimeStep = 0.01;
        assert(numElecTimeStepsPerMechTimestep>0);
        assert(numElementsPerDimInMechanicsMesh>0);
        mNumElementsPerDimInMechanicsMesh = numElementsPerDimInMechanicsMesh;
        mNumElecTimestepsPerMechTimestep = numElecTimeStepsPerMechTimestep;
        mMechanicsTimeStep = mElectricsTimeStep*mNumElecTimestepsPerMechTimestep;
        assert(nhsOdeTimeStep <= mMechanicsTimeStep+1e-14);
        mNhsOdeTimeStep = nhsOdeTimeStep;

        // check whether output is required
        mWriteOutput = (outputDirectory!="");
        if(mWriteOutput)
        {
            mOutputDirectory = outputDirectory;
            mDeformationOutputDirectory = mOutputDirectory + "/deformation";
            HeartConfig::Instance()->SetOutputDirectory(mOutputDirectory + "/electrics");
            HeartConfig::Instance()->SetOutputFilenamePrefix("voltage");
        }
        else
        {
            mDeformationOutputDirectory = ""; // not really necessary but a bit safer, as passed in ConstructMechanicsAssembler
        }
        mNoElectricsOutput = false;

        // initialise all the pointers
        mpElectricsMesh = NULL;
        mpMechanicsMesh = NULL;
        mpCardiacMechAssembler = NULL;

        // Create the Logfile (note we have to do this after the output dir has been
        // created, else the log file might get cleaned away
        std::string log_dir = mOutputDirectory; // just the TESTOUTPUT dir if mOutputDir="";
        LogFile::Instance()->Set(2, mOutputDirectory);
        LogFile::Instance()->WriteHeader("Electromechanics");
        LOG(2, DIM << "d Implicit CardiacElectroMechanics Simulation:");
        LOG(2, "End time = " << mEndTime << ", electrics time step = " << mElectricsTimeStep << ", mechanics timestep = " << mMechanicsTimeStep << "\n");
        LOG(2, "Nhs ode timestep " << mNhsOdeTimeStep);
        LOG(2, "Output is written to " << mOutputDirectory << "/[deformation/electrics]");

        // by default we don't use the direct solver, as not all machines are
        // set up to use UMFPACK yet. However, it is MUCH better than GMRES.
        mUseDirectLinearSolver = false;

        mIsWatchedLocation = false;
        mWatchedElectricsNodeIndex = UNSIGNED_UNSET;
        mWatchedMechanicsNodeIndex = UNSIGNED_UNSET;
    }

    /**
     *  Delete allocated memory and close the watched location file
     */
    virtual ~CardiacElectroMechanicsProblem2()
    {
        delete mpMonodomainProblem;
        delete mpCardiacMechAssembler;
        delete mpElectricsMesh;
        delete mpMechanicsMesh;

        if(mIsWatchedLocation)
        {
            mpWatchedLocationFile->close();
        }

        LogFile::Close();
    }

    /**
     *  Initialise the class. Calls ConstructMeshes() and 
     *  ConstructMechanicsAssembler(). Initialises the MonodomainProblem 
     *  and sets up the electrics mesh to mechanics mesh data.
     */
    void Initialise()
    {
        LOG(2, "Initialising meshes and cardiac mechanics assembler..");

        assert(mpElectricsMesh==NULL);
        assert(mpMechanicsMesh==NULL);
        assert(mpCardiacMechAssembler==NULL);

        // construct the two meshes
        ConstructMeshes();

        if(mIsWatchedLocation)
        {
            DetermineWatchedNodes();
        }

        // initialise monodomain problem
        mpMonodomainProblem->SetMesh(mpElectricsMesh);
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75,1.75,1.75));

        mpMonodomainProblem->Initialise();

        // construct mechanics assembler
        ConstructMechanicsAssembler();

//        if(mUseDirectLinearSolver)
//        {
//            mpCardiacMechAssembler->UseDirectSolver();
//        }

        // find the element nums and weights for each gauss point in the mechanics mesh
        mElementAndWeightsForQuadPoints.resize(mpCardiacMechAssembler->GetTotalNumQuadPoints());

        // get the quad point positions in the mechanics assembler
        QuadraturePointsGroup<DIM> quad_point_posns(*mpMechanicsMesh, *(mpCardiacMechAssembler->GetQuadratureRule()));

        // find the electrics element and weight for each quad point in the mechanics mesh,
        // and store
        for(unsigned i=0; i<quad_point_posns.Size(); i++)
        {
            ChastePoint<DIM> point;

            for(unsigned j=0;j<DIM;j++)
            {
                point.rGetLocation()[j]=quad_point_posns.Get(i)[j];
            }

            unsigned elem_index = mpElectricsMesh->GetContainingElementIndex(point);
            c_vector<double,DIM+1> weight = mpElectricsMesh->GetElement(elem_index)->CalculateInterpolationWeights(point);

            mElementAndWeightsForQuadPoints[i].ElementNum = elem_index;
            mElementAndWeightsForQuadPoints[i].Weights = weight;
        }

        if(mWriteOutput)
        {
            TrianglesMeshWriter<DIM,DIM> mesh_writer(mOutputDirectory,"electrics_mesh",false);
            mesh_writer.WriteFilesUsingMesh(*mpElectricsMesh);
        }

//        // get the assembler to compute which electrics nodes are in each mechanics mesh
//        mpCardiacMechAssembler->ComputeElementsContainingNodes(mpElectricsMesh);
//        assert(DIM==2);
//
//        NodewiseData<DIM>::Instance()->AllocateMemory(mpElectricsMesh->GetNumNodes(), 3);
//        std::vector<std::vector<double> >& r_c_inverse = NodewiseData<DIM>::Instance()->rGetData();
//        for(unsigned i=0; i<r_c_inverse.size(); i++)
//        {
//            r_c_inverse[i][0] = 1.0;
//            r_c_inverse[i][1] = 0.0;
//            r_c_inverse[i][2] = 1.0;
//        }
    }

    /**
     *  Solve the electromechanics problem
     */
    void Solve()
    {
        // initialise the meshes and mechanics assembler
        if(mpCardiacMechAssembler==NULL)
        {
            Initialise();
        }

        BoundaryConditionsContainer<DIM,DIM,1> bcc;
        bcc.DefineZeroNeumannOnMeshBoundary(mpElectricsMesh, 0);
        mpMonodomainProblem->SetBoundaryConditionsContainer(&bcc);

        // get an electrics assembler from the problem. Note that we don't call
        // Solve() on the CardiacProblem class, we do the looping here.
        AbstractDynamicAssemblerMixin<DIM,DIM,1>* p_electrics_assembler
           = mpMonodomainProblem->CreateAssembler();

        // set up initial voltage etc
        Vec voltage=NULL; //This will be set and used later
        Vec initial_voltage = mpMonodomainProblem->CreateInitialCondition();

        unsigned num_quad_points = mpCardiacMechAssembler->GetTotalNumQuadPoints();
        std::vector<NhsCellularMechanicsOdeSystem> cellmech_systems;
        std::vector<double> intracellular_Ca(num_quad_points, 0.0);

        // write the initial position
        unsigned counter = 0;

        TimeStepper stepper(0.0, mEndTime, mMechanicsTimeStep);

        unsigned mech_writer_counter = 0;
        if (mWriteOutput)
        {
            mpCardiacMechAssembler->SetWriteOutput();
            mpCardiacMechAssembler->WriteOutput(mech_writer_counter);

            if(!mNoElectricsOutput)
            {
                mpMonodomainProblem->InitialiseWriter();
                mpMonodomainProblem->WriteOneStep(stepper.GetTime(), initial_voltage);
            }

            if(mIsWatchedLocation)
            {
                WriteWatchedLocationData(stepper.GetTime(), initial_voltage);
            }
        }

        while (!stepper.IsTimeAtEnd())
        {
            LOG(2, "\nCurrent time = " << stepper.GetTime());
            std::cout << "\n\n ** Current time = " << stepper.GetTime();

            LOG(2, "  Solving electrics");

            for(unsigned i=0; i<mNumElecTimestepsPerMechTimestep; i++)
            {
                double current_time = stepper.GetTime() + i*mElectricsTimeStep;
                double next_time = stepper.GetTime() + (i+1)*mElectricsTimeStep;

                // solve the electrics
                p_electrics_assembler->SetTimes(current_time, next_time, mElectricsTimeStep);
                p_electrics_assembler->SetInitialCondition( initial_voltage );

                voltage = p_electrics_assembler->Solve();

                PetscReal min_voltage, max_voltage;
                VecMax(voltage,PETSC_NULL,&max_voltage); //the second param is where the index would be returned
                VecMin(voltage,PETSC_NULL,&min_voltage);
                if(i==0)
                {
                    LOG(2, "  minimum and maximum voltage is " << min_voltage <<", "<<max_voltage);
                }
                else if(i==1)
                {
                    LOG(2, "  ..");
                }

                VecDestroy(initial_voltage);
                initial_voltage = voltage;
            }


//            p_electrics_assembler->SetMatrixIsNotAssembled();

            // compute Ca_I at each quad point (by interpolation, using the info on which
            // electrics element the quad point is in. Then set Ca_I on the mechanics solver
            LOG(2, "  Interpolating Ca_I");
            for(unsigned i=0; i<mElementAndWeightsForQuadPoints.size(); i++)
            {
                double interpolated_Ca_I = 0;

                Element<DIM,DIM>& element = *(mpElectricsMesh->GetElement(mElementAndWeightsForQuadPoints[i].ElementNum));
                for(unsigned node_index = 0; node_index<element.GetNumNodes(); node_index++)
                {
                    unsigned global_node_index = element.GetNodeGlobalIndex(node_index);
                    double Ca_I_at_node = mpMonodomainProblem->GetPde()->GetCardiacCell(global_node_index)->GetIntracellularCalciumConcentration();
                    interpolated_Ca_I += Ca_I_at_node*mElementAndWeightsForQuadPoints[i].Weights(node_index);
                }

                intracellular_Ca[i] = interpolated_Ca_I;
            }

            LOG(2, "  Setting Ca_I. max value = " << Max(intracellular_Ca));

            // NOTE: HERE WE SHOULD REALLY CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
            // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS.. (esp for implicit)
            
            // set [Ca]
            mpCardiacMechAssembler->SetIntracellularCalciumConcentrations(intracellular_Ca);


            // solve the mechanics
            LOG(2, "  Solving mechanics ");
            //double timestep = std::min(0.01, stepper.GetNextTime()-stepper.GetTime());
            mpCardiacMechAssembler->SetWriteOutput(false);
            mpCardiacMechAssembler->Solve(stepper.GetTime(), stepper.GetNextTime(), mNhsOdeTimeStep);

            LOG(2, "    Number of newton iterations = " << mpCardiacMechAssembler->GetNumNewtonIterations());

            // update the current time
            stepper.AdvanceOneTimeStep();
            counter++;

            // output the results
            if(mWriteOutput && (counter%WRITE_EVERY_NTH_TIME==0))
            {
                LOG(2, "  Writing output");
                // write deformed position
                mech_writer_counter++;
                mpCardiacMechAssembler->SetWriteOutput();
                mpCardiacMechAssembler->WriteOutput(mech_writer_counter);

                if(!mNoElectricsOutput)
                {
                    mpMonodomainProblem->mpWriter->AdvanceAlongUnlimitedDimension();
                    mpMonodomainProblem->WriteOneStep(stepper.GetTime(), voltage);
                }

                if(mIsWatchedLocation)
                {
                    WriteWatchedLocationData(stepper.GetTime(), voltage);
                }
            }

//            // setup the Cinverse data;
//            std::vector<std::vector<double> >& r_c_inverse = NodewiseData<DIM>::Instance()->rGetData();
//            mpCardiacMechAssembler->CalculateCinverseAtNodes(mpElectricsMesh, r_c_inverse);
//
//            // write lambda
//            std::stringstream file_name;
//            file_name << "lambda_" << mech_writer_counter << ".dat";
//            mpCardiacMechAssembler->WriteLambda(mOutputDirectory,file_name.str());


            // write the total elapsed time..
            LogFile::Instance()->WriteElapsedTime("  ");
        }

        if ((mWriteOutput) && (!mNoElectricsOutput))
        {
//            if( PetscTools::AmMaster() ) // ie only if master process and results files were written
//            {
//                // call shell script which converts the data to meshalyzer format
//                std::string chaste_2_meshalyzer;
//                std::stringstream space_dim;
//                space_dim << DIM;
//
//                std::string mesh_full_path =   OutputFileHandler::GetChasteTestOutputDirectory()
//                                             + mOutputDirectory + "/electrics_mesh";
//
//                chaste_2_meshalyzer = "anim/chaste2meshalyzer "     // the executable.
//                                      + space_dim.str() + " "       // argument 1 is the dimension.
//                                      + mesh_full_path + " "        // arg 2 is mesh prefix
//                                      + mOutputDirectory + "/electrics/"
//                                      + "voltage "                  // arg 3 is the results folder and prefix,
//                                                                    // relative to the testoutput folder.
//                                      + "last_simulation";          // arg 4 is the output prefix, relative to
//                                                                    // anim folder.
//                system(chaste_2_meshalyzer.c_str());
//            }
//
            HeartConfig::Instance()->Reset();
            mpMonodomainProblem->mpWriter->Close();
            delete mpMonodomainProblem->mpWriter;
        }

        delete p_electrics_assembler;
    }



    // short helper function - the max of a std::vec. this is in the wrong place
    double Max(std::vector<double>& vec)
    {
        double max = -1e200;
        for(unsigned i=0; i<vec.size(); i++)
        {
            if(vec[i]>max) max=vec[i];
        }
        return max;
    }

    /** Call to not write out voltages */
    void SetNoElectricsOutput()
    {
        mNoElectricsOutput = true;
    }

//    /** Use the direct solver when solving linear systems in the
//     *  mechanics. DEFINITELY should be used in experimental work.
//     */
//    void UseDirectLinearSolver()
//    {
//        mUseDirectLinearSolver = true;
//    }

    /**
     *  Set a location to be watched - for which lots of output
     *  is given. Should correspond to nodes in both meshes.
     *
     *  The watched file will have rows that look like:
     *  time x_pos y_pos [z_pos] voltage Ca_i_conc.
     *  
     *  NOTE: for the Calcium - assumes LUO_RUDY IS USED
     */
    void SetWatchedPosition(c_vector<double,DIM> watchedLocation)
    {
        mIsWatchedLocation = true;
        mWatchedLocation = watchedLocation;
    }
};



#endif /*CARDIACELECTROMECHANICSPROBLEM2_HPP_*/
