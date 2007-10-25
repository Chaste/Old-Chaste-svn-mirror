#ifndef ABSTRACTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define ABSTRACTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include "MonodomainProblem.hpp"
#include "AbstractCardiacMechanicsAssembler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "NhsCellularMechanicsOdeSystem.hpp"
#include "FiniteElasticityTools.hpp"
#include "AbstractElasticityAssembler.hpp"
#include "TrianglesMeshWriter.cpp"
#include "LogFile.hpp"

/* todos:
 * 
 * add comments
 * add tests
 * 
 * det F
 * 
 * move mesh stuff out
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
 *  AbstractCardiacElectroMechanicsProblem
 * 
 *  Main class for solved full electro-mechanical problems. Currently subclasses just
 *  define which meshes to use and which assembler to use.
 * 
 *  Solves a monodomain problem (diffusion plus cell models) on a (fine) electrics (chaste) 
 *  mesh,and a mechanics problem (finite elasticity plus NHS cell models) on a coarse (dealii)
 *  mesh. Variable timesteps to implemented soon. An explicit (unstable) or implicit (Jon
 *  Whiteley's algorithm) can be used. 
 * 
 *  The explicit algorithm:
 *  
 *  Store the position in the electrics mesh of each quad point in the mechanics mesh
 *  For every time: 
 *    Solve the monodomain problem (ie integrate ODEs, solve PDE)
 *    Get intracellular [Ca] at each electrics node and interpolate on each mechanics quad point
 *    Set [Ca], current fibre stretch and stretch rate at each mechanics quad point
 *    Integrate NHS models (one for each quad point) explicity
 *    Get active tension at each quad point and set on the mechanics assembler
 *    Solve static finite elasticity (using active tension as a constant 'source')
 *  end
 * 
 *  TODO: alter monodomain equation for the deformation   
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
 */ 
template<unsigned DIM>
class AbstractCardiacElectroMechanicsProblem
{
friend class TestCardiacElectroMechanicsProblem;

protected :
    /*< The cardiac problem class */
    MonodomainProblem<DIM>* mpMonodomainProblem;
    /*< The mechanics assembler */
    AbstractCardiacMechanicsAssembler<DIM>* mpCardiacMechAssembler;  

    /*< End time. The start time is assumed to be 0.0 */
    double mEndTime;
    /*< The timestep. TODO: different timesteps for different bits */
    double mTimeStep;    
    
    /*< A chaste mesh for the electrics */
    ConformingTetrahedralMesh<DIM,DIM>* mpElectricsMesh;
    /*<  A dealii mesh for the mechanics */
    Triangulation<DIM>*                 mpMechanicsMesh;

    /** 
     *  The (electrics-mesh) element numbers saying which element each 
     *  (mechanics-mesh) gauss point is in, and the weight of that gauss point 
     *  for that particular element.
     */
    std::vector<ElementAndWeights<DIM> > mElementAndWeightsForQuadPoints;

    /*< Whether to use an explicit or implicit method */
    bool mUseExplicitMethod;

    /*< Output directory, relative to TEST_OUTPUT */
    std::string mOutputDirectory;
    std::string mDeformationOutputDirectory;
    /*< Whether to write any output */
    bool mWriteOutput;
    /*< when to write output */    
    const static int WRITE_EVERY_NTH_TIME = 1; 
    
    /**
     *  A pure method constructing the mechanics assembler 
     *  @param mechanicsOutputDir The output directory the assembler
     *  should be created with.
     *  */
    virtual void ConstructMechanicsAssembler(std::string mechanicsOutputDir)=0;
    /*< A pure method to be implemented in the concrete class constructing the meshes */
    virtual void ConstructMeshes()=0;
    
    virtual void PostSolve(double currentTime)
    {
    }
    
public :
    /**
     *  Constructor
     *  @param pCellFactory Pointer to a cell factory for the MonodomainProblem class
     *  @param endTime end time, with start time assumed to be 0.
     *  @param timeStep Time step.
     *  @param useExplicitMethod Whether to use an explicit or implicit method
     *  @param outputDirectory. Defaults to "", in which case no output is written
     */
    AbstractCardiacElectroMechanicsProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                           double endTime,
                                           double timeStep,
                                           bool useExplicitMethod,
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
        mTimeStep = timeStep;
        
        // check whether output is required
        mWriteOutput = (outputDirectory!="");
        if(mWriteOutput)
        {
            mOutputDirectory = outputDirectory;
            mDeformationOutputDirectory = mOutputDirectory + "/deformation";
            mpMonodomainProblem->SetOutputDirectory(mOutputDirectory + "/electrics");
            mpMonodomainProblem->SetOutputFilenamePrefix("voltage");
        }
        else
        {
            mDeformationOutputDirectory = ""; // not really necessary but a bit safer, as passed in ConstructMechanicsAssembler
        }
                
        mUseExplicitMethod = useExplicitMethod;
        
        // initialise all the pointers
        mpElectricsMesh = NULL;
        mpMechanicsMesh = NULL;
        mpCardiacMechAssembler = NULL;
        
        // Create the Logfile (note we have to do this after the output dir has been 
        // created, else the log file might get cleaned away
        std::string log_dir = mOutputDirectory; // just the TESTOUTPUT dir if mOutputDir="";
        LogFile::Instance()->Set(1, mOutputDirectory);
        LogFile::Instance()->WriteHeader("Electromechanics");
        LOG(1, DIM << "d CardiacElectroMechanics Simulation:");
        LOG(1, "End time = " << mEndTime << ", timestep = " << mTimeStep <<  "\n");
        LOG(1, "Output is written to " << mOutputDirectory << "/[deformation/electrics]");
        
        if(mUseExplicitMethod)
        {
            LOG(1, "Solving with explicit method..");
        }
        else
        {
            LOG(1, "Solving with implicit method..");
        }        
    }   
    
    virtual ~AbstractCardiacElectroMechanicsProblem()
    {
        delete mpMonodomainProblem;
        delete mpCardiacMechAssembler;
        delete mpElectricsMesh;
        delete mpMechanicsMesh;
        
        LogFile::Close();
    }
    
    
    /**
     *  Initialise the class. Calls ConstructMeshes() and ConstructMechanicsAssembler() on
     *  the concrete classes. Initialises the MonodomainProblem and sets up the electrics 
     *  mesh to mechanics mesh data.
     */
    void Initialise()
    {
        LOG(1, "Initialising meshes and cardiac mechanics assembler..");
        
        assert(mpElectricsMesh==NULL);
        assert(mpMechanicsMesh==NULL);
        assert(mpCardiacMechAssembler==NULL);
        
        // construct the two meshes
        ConstructMeshes();     

        // initialise monodomain problem                        
        mpMonodomainProblem->SetMesh(mpElectricsMesh);
        mpMonodomainProblem->Initialise();

        // construct mechanics assembler 
        ConstructMechanicsAssembler(mDeformationOutputDirectory);

        // find the element nums and weights for each gauss point in the mechanics mesh
        mElementAndWeightsForQuadPoints.resize(mpCardiacMechAssembler->GetTotalNumQuadPoints());

        // get the quad point positions in the mechanics assembler
        std::vector<std::vector<double> > quad_point_posns
           = FiniteElasticityTools<DIM>::GetQuadPointPositions(*mpMechanicsMesh, mpCardiacMechAssembler->GetNumQuadPointsInEachDimension());

        // find the electrics element and weight for each quad point in the mechanics mesh,
        // and store
        for(unsigned i=0; i<quad_point_posns.size(); i++)
        {
            ChastePoint<DIM> point;

            for(unsigned j=0;j<DIM;j++)
            {
                point.rGetLocation()[j]=quad_point_posns[i][j];
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
    }

    /** 
     *  Solve the electromechanincs problem
     */    
    void Solve()
    {
        // initialise the meshes and mechanics assembler
        if(mpCardiacMechAssembler==NULL)
        {
            Initialise();
        }
        
        // get an electrics assembler from the problem. Note that we don't call
        // Solve() on the CardiacProblem class, we do the looping here.
        AbstractDynamicAssemblerMixin<DIM,DIM,1>* p_electrics_assembler 
           = mpMonodomainProblem->CreateAssembler();

        // set up initial voltage etc
        Vec voltage;        
        Vec initial_voltage = mpMonodomainProblem->CreateInitialCondition();

        // Create stores of lambda, lambda_dot and old lambda
        // Note: these are only needed if an explicit method is used
        unsigned num_quad_points = mpCardiacMechAssembler->GetTotalNumQuadPoints();
        std::vector<double> lambda;
        std::vector<double> old_lambda;
        std::vector<double> dlambda_dt;
        std::vector<NhsCellularMechanicsOdeSystem> cellmech_systems;
        EulerIvpOdeSolver euler_solver;

        // this is the active tension if explicit and the calcium conc if implicit
        std::vector<double> forcing_quantity(num_quad_points, 0.0);
        
        // initial cellmechanics systems, lambda, etc, if required
        if(mUseExplicitMethod)
        {
            lambda.resize(num_quad_points, 1.0);
            old_lambda.resize(num_quad_points, 1.0);
            dlambda_dt.resize(num_quad_points, 0.0);
            cellmech_systems.resize(num_quad_points);
        }
            

        // write the initial position
        // NOTE: small architecture issue here. All concrete assembler (at the moment) are
        // AbstractElasticityAssembler assemblers (naturally) as well as AbstractCardiacMech
        // assemblers, but the compiler doesn't know that here, hence the cast
        
        unsigned counter = 0;

        TimeStepper stepper(0.0, mEndTime, mTimeStep);

        unsigned mech_writer_counter = 0;
        if (mWriteOutput)
        {
            dynamic_cast<AbstractElasticityAssembler<DIM>*>(mpCardiacMechAssembler)->WriteOutput(mech_writer_counter);
            mpMonodomainProblem->InitialiseWriter();
            mpMonodomainProblem->WriteOneStep(stepper.GetTime(), initial_voltage);
        }

        while (!stepper.IsTimeAtEnd())
        {
            LOG(1, "\nCurrent time = " << stepper.GetTime());
            
            // solve the electrics
            p_electrics_assembler->SetTimes(stepper.GetTime(), stepper.GetNextTime(), mTimeStep);
            p_electrics_assembler->SetInitialCondition( initial_voltage );
            
            LOG(1, "  Solving electrics");
            voltage = p_electrics_assembler->Solve();
            
            PetscReal min_voltage, max_voltage;
            VecMax(voltage,PETSC_NULL,&max_voltage); //the second param is where the index would be returned
            VecMin(voltage,PETSC_NULL,&min_voltage);
            LOG(1, "  minimun and maximum voltage is " << min_voltage <<", "<<max_voltage);

            VecDestroy(initial_voltage);
            initial_voltage = voltage;
            
            // compute Ca_I at each quad point (by interpolation, using the info on which
            // electrics element the quad point is in. Then: 
            //   Explicit: Set Ca_I on the nhs systems and solve them to get the active tension
            //   Implicit: Set Ca_I on the mechanics solver

            LOG(1, "  Interpolating Ca_I\n  (and solving NHS models if explicit)");
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

                if(mUseExplicitMethod)
                {
                    // explicit: forcing quantity on the assembler is the active tension
                    cellmech_systems[i].SetLambdaAndDerivative(lambda[i], dlambda_dt[i]);
                    cellmech_systems[i].SetIntracellularCalciumConcentration(interpolated_Ca_I);
        
                    euler_solver.SolveAndUpdateStateVariable(&cellmech_systems[i], stepper.GetTime(), stepper.GetNextTime(), mTimeStep);
                    forcing_quantity[i] = cellmech_systems[i].GetActiveTension();
                }
                else
                {
                    // explicit: forcing quantity on the assembler is the calcium concentration
                    forcing_quantity[i] = interpolated_Ca_I;
                }
            }

            if(mUseExplicitMethod)
            {
                LOG(1, "  Setting active tension. max value = " << Max(forcing_quantity));
            }
            else
            {
                LOG(1, "  Setting Ca_I. max value = " << Max(forcing_quantity));
            }

            // NOTE: HERE WE SHOULD REALLY CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
            // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS.. (esp for implicit)
            
            // set the active tensions if explicit, or the [Ca] if implicit
            mpCardiacMechAssembler->SetForcingQuantity(forcing_quantity);

            // solve the mechanics
            LOG(1, "  Solving mechanics");
            mpCardiacMechAssembler->Solve(stepper.GetTime(), stepper.GetNextTime(), stepper.GetNextTime()-stepper.GetTime());
            
            unsigned num_iters = dynamic_cast<AbstractElasticityAssembler<DIM>*>(mpCardiacMechAssembler)->GetNumNewtonIterations();
            LOG(1, "    Number of newton iterations = " << num_iters);

            // if explicit store the new lambda and update lam
            if(mUseExplicitMethod)
            {
                // update lambda and dlambda_dt;
                LOG(1, "  Updating lambda");
                old_lambda = lambda;
                lambda = mpCardiacMechAssembler->rGetLambda();
                for(unsigned i=0; i<dlambda_dt.size(); i++)
                {
                    dlambda_dt[i] = (lambda[i] - old_lambda[i])/mTimeStep;
                }
            }

            if(mWriteOutput && (counter%WRITE_EVERY_NTH_TIME==0))
            {
                LOG(1, "  Writing output");
                // write deformed position                    
                mech_writer_counter++;
                dynamic_cast<AbstractElasticityAssembler<DIM>*>(mpCardiacMechAssembler)->WriteOutput(mech_writer_counter);

                mpMonodomainProblem->mpWriter->AdvanceAlongUnlimitedDimension();
                mpMonodomainProblem->WriteOneStep(stepper.GetTime(), voltage);
            }
            
            PostSolve(stepper.GetTime());
            
            //// TODO: Update the monodomain equations for the deformation
            
            // update the current time
            stepper.AdvanceOneTimeStep();
            counter++;

            // write the total elapsed time..
            LogFile::Instance()->WriteElapsedTime("  ");
        }


        if (mWriteOutput)
        {
            if(mpMonodomainProblem->mpWriter->AmMaster()) // ie only if master process and results files were written
            {
                // call shell script which converts the data to meshalyzer format
                std::string chaste_2_meshalyzer;
                std::stringstream space_dim;
                space_dim << DIM;
                
                std::string mesh_full_path =   OutputFileHandler::GetChasteTestOutputDirectory()
                                             + mOutputDirectory + "/electrics_mesh";
                
                chaste_2_meshalyzer = "anim/chaste2meshalyzer "     // the executable.
                                      + space_dim.str() + " "       // argument 1 is the dimension.
                                      + mesh_full_path + " "        // arg 2 is mesh prefix
                                      + mOutputDirectory + "/electrics/"
                                      + "voltage "                  // arg 3 is the results folder and prefix,
                                                                    // relative to the testoutput folder.
                                      + "last_simulation";          // arg 4 is the output prefix, relative to
                                                                    // anim folder.                
                system(chaste_2_meshalyzer.c_str());
            }

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
};



#endif /*ABSTRACTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
