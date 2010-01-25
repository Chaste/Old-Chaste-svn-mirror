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

#include "CardiacElectroMechanicsProblem.hpp"

#include "OutputFileHandler.hpp"
#include "ReplicatableVector.hpp"
#include "HeartConfig.hpp"
#include "LogFile.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "TimeStepper.hpp"
#include "QuadraturePointsGroup.hpp"
#include "TrianglesMeshWriter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "PetscTools.hpp"

#include "ImplicitCardiacMechanicsAssembler.hpp"
#include "ExplicitCardiacMechanicsAssembler.hpp"
// if including Cinv in monobidomain equations
//#include "NodewiseData.hpp"

#include "MooneyRivlinMaterialLaw.hpp"

#include "CmguiDeformedSolutionsWriter.hpp"


template<unsigned DIM>
void CardiacElectroMechanicsProblem<DIM>::DetermineWatchedNodes()
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
        #define COVERAGE_IGNORE
        std::cout << "ERROR: Could not find an electrics node very close to requested watched location - "
                  << "min distance was " << min_dist << " for node " << node_index
                  << " at location " << pos << std::flush;;

        //// the following causes a seg fault for some reason (!!???!!!)
        //EXCEPTION("Could not find an electrics node very close to requested watched location");
        NEVER_REACHED;
        #undef COVERAGE_IGNORE
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
        #define COVERAGE_IGNORE
        std::cout << "ERROR: Could not find a mechanics node very close to requested watched location - "
                  << "min distance was " << min_dist << " for node " << node_index
                  << " at location " << pos_at_min;

        //// the following causes a seg fault for some reason (!!???!!!)
        //EXCEPTION("Could not find a mechanics node very close to requested watched location");
        assert(0);
        #undef COVERAGE_IGNORE
    }
    else
    {
        LOG_AND_COUT(2,"Chosen electrics node "<<node_index<<" at location " << pos << " to be watched");
        mWatchedMechanicsNodeIndex = node_index;
    }

    OutputFileHandler handler(mOutputDirectory);
    mpWatchedLocationFile = handler.OpenOutputFile("watched.txt");
}

template<unsigned DIM>
void CardiacElectroMechanicsProblem<DIM>::WriteWatchedLocationData(double time, Vec voltage)
{
    assert(mIsWatchedLocation);

    std::vector<c_vector<double,DIM> >& deformed_position = mpCardiacMechAssembler->rGetDeformedPosition();

    ///\todo Improve efficiency of this method?
    ReplicatableVector voltage_replicated(voltage);
    double V = voltage_replicated[mWatchedElectricsNodeIndex];

    //// Removed the following which also took this nodes calcium concentration and printing, because (it isn't that 
    //// important and) it won't work in parallel and has the hardcoded index issue described below. 
    //     // \todo: NOTE!!! HARDCODED state variable index - assumes Lr91. 
    //     // Metadata is currently being added to CellML models and then this will be avoided by asking for Calcium.
    //    double Ca = mpMonodomainProblem->GetMonodomainPde()->GetCardiacCell(mWatchedElectricsNodeIndex)->rGetStateVariables()[3];

    *mpWatchedLocationFile << time << " ";
    for(unsigned i=0; i<DIM; i++)
    {
        *mpWatchedLocationFile << deformed_position[mWatchedMechanicsNodeIndex](i) << " ";
    }
    *mpWatchedLocationFile << V << "\n";
    mpWatchedLocationFile->flush();
}


template<unsigned DIM>
CardiacElectroMechanicsProblem<DIM>::CardiacElectroMechanicsProblem(
            ContractionModel contractionModel,
            TetrahedralMesh<DIM,DIM>* pElectricsMesh,
            QuadraticMesh<DIM>* pMechanicsMesh,
            std::vector<unsigned> fixedMechanicsNodes,
            AbstractCardiacCellFactory<DIM>* pCellFactory,
            double endTime,
            double electricsPdeTimeStep,
            unsigned numElecTimeStepsPerMechTimestep,
            double contractionModelOdeTimeStep,
            std::string outputDirectory = "")
{
    // Start-up mechanics event handler..
    MechanicsEventHandler::Reset();
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ALL);
    // disable the electric event handler, because we use a problem class but
    // don't call Solve, so we would have to worry about starting and ending any
    // events in AbstractCardiacProblem::Solve() (esp. calling EndEvent(EVERYTHING))
    // if we didn't disable it.
    HeartEventHandler::Disable();
    
    mContractionModel = contractionModel;

    // create the monodomain problem. Note the we use this to set up the cells,
    // get an initial condition (voltage) vector, and get an assembler. We won't
    // ever call solve on the MonodomainProblem
    assert(pCellFactory != NULL);
    mpMonodomainProblem = new MonodomainProblem<DIM>(pCellFactory);

    // save time infomation
    assert(endTime > 0);
    mEndTime = endTime;

    assert(electricsPdeTimeStep>0);
    mElectricsTimeStep = electricsPdeTimeStep;

    assert(numElecTimeStepsPerMechTimestep>0);

    mNumElecTimestepsPerMechTimestep = numElecTimeStepsPerMechTimestep;

    mMechanicsTimeStep = mElectricsTimeStep*mNumElecTimestepsPerMechTimestep;
    assert(contractionModelOdeTimeStep <= mMechanicsTimeStep+1e-14);
    mContractionModelOdeTimeStep = contractionModelOdeTimeStep;

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
        mDeformationOutputDirectory = "";
    }
    mNoElectricsOutput = false;

    // initialise all the pointers
    mpElectricsMesh = pElectricsMesh; // note these are allowed to be null, in case a child constructor wants to create them
    mpMechanicsMesh = pMechanicsMesh;
    mFixedNodes = fixedMechanicsNodes;

    mpCardiacMechAssembler = NULL;

    // Create the Logfile (note we have to do this after the output dir has been
    // created, else the log file might get cleaned away
    std::string log_dir = mOutputDirectory; // just the TESTOUTPUT dir if mOutputDir="";
    LogFile::Instance()->Set(2, mOutputDirectory);
    LogFile::Instance()->WriteHeader("Electromechanics");
    LOG(2, DIM << "d Implicit CardiacElectroMechanics Simulation:");
    LOG(2, "End time = " << mEndTime << ", electrics time step = " << mElectricsTimeStep << ", mechanics timestep = " << mMechanicsTimeStep << "\n");
    LOG(2, "Contraction model ode timestep " << mContractionModelOdeTimeStep);
    LOG(2, "Output is written to " << mOutputDirectory << "/[deformation/electrics]");
#define COVERAGE_IGNORE
/// \todo Cover these lines
    if(mpElectricsMesh != NULL)
    {
        LOG(2, "Electrics mesh has " << mpElectricsMesh->GetNumNodes() << " nodes");
    }
    if(mpMechanicsMesh != NULL)
    {
        LOG(2, "Mechanics mesh has " << mpMechanicsMesh->GetNumNodes() << " nodes");
    }
#undef COVERAGE_IGNORE
    mIsWatchedLocation = false;
    mWatchedElectricsNodeIndex = UNSIGNED_UNSET;
    mWatchedMechanicsNodeIndex = UNSIGNED_UNSET;
    
    mFibreDirectionsFile = "";
}

template<unsigned DIM>
CardiacElectroMechanicsProblem<DIM>::~CardiacElectroMechanicsProblem()
{   /** 
     * NOTE if SetWatchedLocation but not Initialise has been
     * called, mpWatchedLocationFile will be uninitialised and
     * using it will cause a seg fault. Hence the mpMechanicsMesh!=NULL
     * it is true if Initialise has been called.
     */
    if(mIsWatchedLocation && mpMechanicsMesh)
    {
        mpWatchedLocationFile->close();
    }

    delete mpMonodomainProblem;

    delete mpCardiacMechAssembler;

    LogFile::Close();
}

template<unsigned DIM>
void CardiacElectroMechanicsProblem<DIM>::Initialise()
{
    LOG(2, "Initialising meshes and cardiac mechanics assembler..");

    assert(mpElectricsMesh!=NULL);
    assert(mpMechanicsMesh!=NULL);
    assert(mpCardiacMechAssembler==NULL);

    if(mIsWatchedLocation)
    {
        DetermineWatchedNodes();
    }


    // initialise monodomain problem
    mpMonodomainProblem->SetMesh(mpElectricsMesh);
    HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75,1.75,1.75));
    mpMonodomainProblem->Initialise();

    // Construct mechanics assembler
    // Here we pick the best solver for each particular contraction model. Commented out versions are 
    // for experimentation.
    switch(mContractionModel)
    {
        case NASH2004:
            // stretch and stretch-rate independent, so should use explicit
            mpCardiacMechAssembler = new ExplicitCardiacMechanicsAssembler<DIM>(mContractionModel,mpMechanicsMesh,mDeformationOutputDirectory,mFixedNodes);
            break;
        case KERCHOFFS2003:
            // stretch independent, so should use implicit (explicit may be unstable)
            mpCardiacMechAssembler = new ImplicitCardiacMechanicsAssembler<DIM>(mContractionModel,mpMechanicsMesh,mDeformationOutputDirectory,mFixedNodes);
            break;
        case NHS:
            // stretch and stretch-rate independent, so should definitely use implicit
            mpCardiacMechAssembler = new ImplicitCardiacMechanicsAssembler<DIM>(mContractionModel,mpMechanicsMesh,mDeformationOutputDirectory,mFixedNodes);
            break;
        default:
            EXCEPTION("Invalid contraction model, options are: KERCHOFFS2003 or NHS");
    }
    
    if(mFibreDirectionsFile!="")
    {
       mpCardiacMechAssembler->SetVariableFibreDirections(mFibreDirectionsFile);
    }
    
    // find the element nums and weights for each gauss point in the mechanics mesh
    mElementAndWeightsForQuadPoints.resize(mpCardiacMechAssembler->GetTotalNumQuadPoints());

    // get the quad point positions in the mechanics assembler
    QuadraturePointsGroup<DIM> quad_point_posns(*mpMechanicsMesh, *(mpCardiacMechAssembler->GetQuadratureRule()));

    // find the electrics element and weight for each quad point in the mechanics mesh,
    // and store
    unsigned last_element = 0;
    for(unsigned i=0; i<quad_point_posns.Size(); i++)
    {
    	//// this can be slow with certains pairs of mesh? #1198
		//std::cout << "\r " << i << " of " << quad_point_posns.Size();
        
        ChastePoint<DIM> point;

        for(unsigned j=0; j<DIM; j++)
        {
            point.rGetLocation()[j]=quad_point_posns.Get(i)[j];
        }

        unsigned elem_index = mpElectricsMesh->GetContainingElementIndexWithInitialGuess(point, last_element);
        last_element = elem_index;
        
        c_vector<double,DIM+1> weight = mpElectricsMesh->GetElement(elem_index)->CalculateInterpolationWeights(point);

        mElementAndWeightsForQuadPoints[i].ElementNum = elem_index;
        mElementAndWeightsForQuadPoints[i].Weights = weight;
    }

    if(mWriteOutput)
    {
        TrianglesMeshWriter<DIM,DIM> mesh_writer(mOutputDirectory,"electrics_mesh",false);
        mesh_writer.WriteFilesUsingMesh(*mpElectricsMesh);
    }

//// when using *Cinverse* in electrics
//    // get the assembler to compute which electrics nodes are in each mechanics mesh
//    mpCardiacMechAssembler->ComputeElementsContainingNodes(mpElectricsMesh);
//    assert(DIM==2);
//
//    NodewiseData<DIM>::Instance()->AllocateMemory(mpElectricsMesh->GetNumNodes(), 3);
//    std::vector<std::vector<double> >& r_c_inverse = NodewiseData<DIM>::Instance()->rGetData();
//    for(unsigned i=0; i<r_c_inverse.size(); i++)
//    {
//        r_c_inverse[i][0] = 1.0;
//        r_c_inverse[i][1] = 0.0;
//        r_c_inverse[i][2] = 1.0;
//    }
}

template<unsigned DIM>
void CardiacElectroMechanicsProblem<DIM>::Solve()
{
    
    // initialise the meshes and mechanics assembler
    if(mpCardiacMechAssembler==NULL)
    {
        Initialise();
    }

    boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>);
    p_bcc->DefineZeroNeumannOnMeshBoundary(mpElectricsMesh, 0);
    mpMonodomainProblem->SetBoundaryConditionsContainer(p_bcc);

    // get an electrics assembler from the problem. Note that we don't call
    // Solve() on the CardiacProblem class, we do the looping here.
    AbstractDynamicAssemblerMixin<DIM,DIM,1>* p_electrics_assembler
       = mpMonodomainProblem->CreateAssembler();

    // set up initial voltage etc
    Vec voltage=NULL; //This will be set and used later
    Vec initial_voltage = mpMonodomainProblem->CreateInitialCondition();

    unsigned num_quad_points = mpCardiacMechAssembler->GetTotalNumQuadPoints();
    std::vector<double> interpolated_calcium_concs(num_quad_points, 0.0);
    std::vector<double> interpolated_voltages(num_quad_points, 0.0);

    // write the initial position
    unsigned counter = 0;

    TimeStepper stepper(0.0, mEndTime, mMechanicsTimeStep);

    CmguiDeformedSolutionsWriter<DIM>* p_cmgui_writer = NULL;

    unsigned mech_writer_counter = 0;
    if (mWriteOutput)
    {
        mpCardiacMechAssembler->SetWriteOutput();
        mpCardiacMechAssembler->WriteOutput(mech_writer_counter);

        p_cmgui_writer = new CmguiDeformedSolutionsWriter<DIM>(mOutputDirectory+"/cmgui", "solution", *(this->mpMechanicsMesh));
        p_cmgui_writer->WriteInitialMesh();


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
        std::cout << "\n\n ** Current time = " << stepper.GetTime() << "\n";

        LOG(2, "  Solving electrics");
        MechanicsEventHandler::BeginEvent(MechanicsEventHandler::NON_MECH);
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

//// when using *Cinverse* in electrics
//        p_electrics_assembler->SetMatrixIsNotAssembled();

        // compute Ca_I at each quad point (by interpolation, using the info on which
        // electrics element the quad point is in. Then set Ca_I on the mechanics solver
        LOG(2, "  Interpolating Ca_I and voltage");
        
        ReplicatableVector voltage_repl(voltage);
        
        // collect all the calcium concentrations (from the cells, which are
        // distributed) in one (replicated) vector
        ReplicatableVector calcium_repl(mpElectricsMesh->GetNumNodes());
        PetscInt lo, hi;
        VecGetOwnershipRange(voltage, &lo, &hi);        
        for(int index=lo; index<hi; index++)
        {
            calcium_repl[index] = mpMonodomainProblem->GetPde()->GetCardiacCell(index)->GetIntracellularCalciumConcentration();
        }
        PetscTools::Barrier();
        calcium_repl.Replicate(lo,hi);
        
        
        for(unsigned i=0; i<mElementAndWeightsForQuadPoints.size(); i++)
        {
            double interpolated_CaI = 0;
            double interpolated_voltage = 0;

            Element<DIM,DIM>& element = *(mpElectricsMesh->GetElement(mElementAndWeightsForQuadPoints[i].ElementNum));
            for(unsigned node_index = 0; node_index<element.GetNumNodes(); node_index++)
            {
                unsigned global_node_index = element.GetNodeGlobalIndex(node_index);
                double CaI_at_node =  calcium_repl[global_node_index]; 
                interpolated_CaI += CaI_at_node*mElementAndWeightsForQuadPoints[i].Weights(node_index);
                interpolated_voltage += voltage_repl[global_node_index]*mElementAndWeightsForQuadPoints[i].Weights(node_index);
            }

            interpolated_calcium_concs[i] = interpolated_CaI;
            interpolated_voltages[i] = interpolated_voltage;
        }

        LOG(2, "  Setting Ca_I. max value = " << Max(interpolated_calcium_concs));

        // NOTE IF NHS: HERE WE SHOULD PERHAPS CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
        // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS.. (esp for implicit)

        // set [Ca], V, t
        mpCardiacMechAssembler->SetCalciumAndVoltage(interpolated_calcium_concs, interpolated_voltages);
        MechanicsEventHandler::EndEvent(MechanicsEventHandler::NON_MECH);

        // solve the mechanics
        LOG(2, "  Solving mechanics ");
        //double timestep = std::min(0.01, stepper.GetNextTime()-stepper.GetTime());
        mpCardiacMechAssembler->SetWriteOutput(false);

        MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ALL_MECH);
        mpCardiacMechAssembler->Solve(stepper.GetTime(), stepper.GetNextTime(), mContractionModelOdeTimeStep);
        MechanicsEventHandler::EndEvent(MechanicsEventHandler::ALL_MECH);

        LOG(2, "    Number of newton iterations = " << mpCardiacMechAssembler->GetNumNewtonIterations());

        // update the current time
        stepper.AdvanceOneTimeStep();
        counter++;

        // output the results
        MechanicsEventHandler::BeginEvent(MechanicsEventHandler::OUTPUT);
        if(mWriteOutput && (counter%WRITE_EVERY_NTH_TIME==0))
        {
            LOG(2, "  Writing output");
            // write deformed position
            mech_writer_counter++;
            mpCardiacMechAssembler->SetWriteOutput();
            mpCardiacMechAssembler->WriteOutput(mech_writer_counter);

            p_cmgui_writer->WriteDeformationPositions(rGetDeformedPosition(), counter);

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
        MechanicsEventHandler::EndEvent(MechanicsEventHandler::OUTPUT);

//// when using *Cinverse* in electrics
//        // setup the Cinverse data;
//        std::vector<std::vector<double> >& r_c_inverse = NodewiseData<DIM>::Instance()->rGetData();
//        mpCardiacMechAssembler->CalculateCinverseAtNodes(mpElectricsMesh, r_c_inverse);
//
//        // write lambda
//        std::stringstream file_name;
//        file_name << "lambda_" << mech_writer_counter << ".dat";
//        mpCardiacMechAssembler->WriteLambda(mOutputDirectory,file_name.str());


        // write the total elapsed time..
        LogFile::Instance()->WriteElapsedTime("  ");
    }

    if ((mWriteOutput) && (!mNoElectricsOutput))
    {
        HeartConfig::Instance()->Reset();
        mpMonodomainProblem->mpWriter->Close();
        delete mpMonodomainProblem->mpWriter;


        // Convert simulation data to Meshalyzer format
        //
        std::string input_dir = mOutputDirectory+"/electrics";
        std::string config_directory = HeartConfig::Instance()->GetOutputDirectory();
        HeartConfig::Instance()->SetOutputDirectory(input_dir);
        Hdf5ToMeshalyzerConverter<DIM,DIM> converter(input_dir, "voltage", mpElectricsMesh);
        
        // Write mesh in a suitable form for meshalyzer
        std::string output_directory =  mOutputDirectory + "/electrics/output";
        // Write the mesh
        MeshalyzerMeshWriter<DIM,DIM> mesh_writer(output_directory, "mesh", false);
        mesh_writer.WriteFilesUsingMesh(*mpElectricsMesh);
        // Write the parameters out
        HeartConfig::Instance()->Write();
        
        // reset to the default value
        HeartConfig::Instance()->SetOutputDirectory(config_directory);
        
    }

    if(p_cmgui_writer)
    {
        p_cmgui_writer->WriteCmguiScript();
        delete p_cmgui_writer;
    }
    VecDestroy(voltage);
    delete p_electrics_assembler;

    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ALL);
}



template<unsigned DIM>
double CardiacElectroMechanicsProblem<DIM>::Max(std::vector<double>& vec)
{
    double max = -1e200;
    for(unsigned i=0; i<vec.size(); i++)
    {
        if(vec[i]>max) max=vec[i];
    }
    return max;
}

template<unsigned DIM>
void CardiacElectroMechanicsProblem<DIM>::SetNoElectricsOutput()
{
    mNoElectricsOutput = true;
}

template<unsigned DIM>
void CardiacElectroMechanicsProblem<DIM>::SetWatchedPosition(c_vector<double,DIM> watchedLocation)
{
    mIsWatchedLocation = true;
    mWatchedLocation = watchedLocation;
}

template<unsigned DIM>
void CardiacElectroMechanicsProblem<DIM>::SetVariableFibreDirectionsFile(std::string fibreDirectionsFile)
{
    mFibreDirectionsFile = fibreDirectionsFile;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& CardiacElectroMechanicsProblem<DIM>::rGetDeformedPosition()
{
    return mpCardiacMechAssembler->rGetDeformedPosition();
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

//note: 1d incompressible material doesn't make sense
template class CardiacElectroMechanicsProblem<2>;
template class CardiacElectroMechanicsProblem<3>;
