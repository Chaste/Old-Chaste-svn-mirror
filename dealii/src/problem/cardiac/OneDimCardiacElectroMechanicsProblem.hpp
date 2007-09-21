#ifndef ONDEDIMCARDIACELECTROMECHANICSPROBLEM_HPP_
#define ONDEDIMCARDIACELECTROMECHANICSPROBLEM_HPP_

#include "MonodomainProblem.hpp"
#include "OneDimCardiacMechanicsAssembler.hpp"
#include "ImplicitOneDimCardiacMechanicsAssembler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "FiniteElasticityTools.hpp"
#include "NhsCellularMechanicsOdeSystem.hpp"


/* todos:
 * 
 * add comments
 * add tests
 * 
 * move mesh stuff out
 * 
 * refactor to use 2d/3d cardiacmechassembler
 * 
 * think about architecture (of AbstractCardiacProblem) when this is done properly..
 */


///\todo -make this struct
/**
 *  At the beginning of a two mesh simulation we need to figure out and store
 *  which (electrics-mesh) element each (mechanics-mesh) gauss point is in, and
 *  what the weight of that gauss point for that particular element is. This struct
 *  just contains this two pieces of data
 */
template<unsigned DIM>
class ElementAndWeights
{
public:
    unsigned ElementNum;
    c_vector<double,DIM+1> Weights;  
};



/**
 *  OneDimCardiacElectroMechanicsProblem 
 *  
 *  Solve a coupled cardiac electromechanics problem
 *  
 *  Currently just uses the monodomain assembler (to use bidomain need to manually change 
 *  the code), and currently just uses the OneDimCardiacMechanicsAssembler, will eventually
 *  use the (2d/3d) CardiacMechanicsAssembler. Can use both explicit and implicit 1d 
 *  assemblers
 */
template<unsigned DIM>
class OneDimCardiacElectroMechanicsProblem
{
private: 
    /*< The cardiac problem class */
    MonodomainProblem<DIM>* mpMonodomainProblem;
    
    /*< The mechanics assembler */
    OneDimCardiacMechanicsAssembler* mpCardiacMechAssembler;  // can be ImplicitOneDim..

    /*< End time. The start time is assumed to be 0.0 */
    double mEndTime;
    /*< The timestep. TODO: different timesteps for different bits */
    double mTimeStep;    
    
    /*< A chaste mesh for the electrics */
    ConformingTetrahedralMesh<DIM,DIM>* mpElectricsMesh;
    /*<  A dealii mesh for the mechanics */
    Triangulation<DIM>*                 mpMechanicsMesh;

    /** 
     *  The (electrics-mesh) element numbers and saying which element each 
     *  (mechanics-mesh) gauss point is in, and the weight of that gauss point 
     *  for that particular element is.
     */
    std::vector<ElementAndWeights<DIM> > mElementAndWeightsForQuadPoints;

    /*< Whether to use an explicit or implicit method */
    bool mUseExplicitMethod;

    /*< Output directory, relative to TEST_OUTPUT */
    std::string mOutputDirectory;
    /*< Whether to write any output */
    bool mWriteOutput;
    /*< when to write output */    
    const static int WRITE_EVERY_NTH_TIME = 1; 
        
public:
    /** 
     *  Constructor
     *  @param pCellFactory cell factory for creating cells (see Monodomain tests)
     *  @endTime end time of the simulation. Start time is assumed to be 0.0
     *  @timeStep time step for the electrics (and currently the mechanics too)
     *  @useExplicit Whether to use an explicit or implicit mechanics solver
     *  @outputDirectory. Output directory. Omit if no output is required.
     * 
     *  The meshes are currently hardcoded in here.
     */
    OneDimCardiacElectroMechanicsProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   double endTime,
                                   double timeStep,
                                   bool useExplicitMethod,
                                   std::string outputDirectory = "")
    {
        assert(pCellFactory != NULL);
        mpMonodomainProblem = new MonodomainProblem<DIM>(pCellFactory);
        
        assert(endTime > 0);
        mEndTime = endTime;
        mTimeStep = timeStep;
        
        // check whether output is required
        mWriteOutput = (outputDirectory!="");
        mOutputDirectory = outputDirectory;
                
        mUseExplicitMethod = useExplicitMethod;        
                        
        // create electrics mesh
        mpElectricsMesh = new ConformingTetrahedralMesh<DIM,DIM>();
        unsigned num_elem = 128;
        mpElectricsMesh->ConstructLinearMesh(num_elem);
        mpElectricsMesh->Scale(1.0/num_elem);
        
        assert(DIM==1);
        mpMonodomainProblem->SetMesh(mpElectricsMesh);
        mpMonodomainProblem->Initialise();
        
        // create mechanics mesh
        mpMechanicsMesh = new Triangulation<DIM>();
        GridGenerator::hyper_cube(*mpMechanicsMesh, 0.0, 1);
        mpMechanicsMesh->refine_global(7);
        
        assert(mpMechanicsMesh->n_vertices()==mpElectricsMesh->GetNumNodes());

        if(mUseExplicitMethod)
        {
            mpCardiacMechAssembler = new OneDimCardiacMechanicsAssembler(mpMechanicsMesh);
        }
        else
        {
            mpCardiacMechAssembler = new ImplicitOneDimCardiacMechanicsAssembler(mpMechanicsMesh);
        }
        
        // find the element nums and weights for each gauss point in the mechanics mesh
        mElementAndWeightsForQuadPoints.resize(mpCardiacMechAssembler->GetTotalNumQuadPoints());
        std::vector<std::vector<double> > quad_point_posns
           = FiniteElasticityTools<DIM>::GetQuadPointPositions(*mpMechanicsMesh, mpCardiacMechAssembler->GetNumQuadPointsInEachDimension());
        
        for(unsigned i=0; i<quad_point_posns.size(); i++)
        {
            assert(DIM==1);
            ChastePoint<DIM> point(quad_point_posns[i][0]);
            
            unsigned elem_index = mpElectricsMesh->GetContainingElementIndex(point);
            c_vector<double,DIM+1> weight = mpElectricsMesh->GetElement(elem_index)->CalculateInterpolationWeights(point);
            
            mElementAndWeightsForQuadPoints[i].ElementNum = elem_index;
            mElementAndWeightsForQuadPoints[i].Weights = weight;
        }
    }
    

    /** 
     *  Solve the electromechanincs problem
     */    
    void Solve()
    {
        // get an electrics assembler from the problem. Note that we don't call
        // Solve() on the CardiacProblem class, we do the looping here.
        AbstractDynamicAssemblerMixin<DIM,DIM,1>* mpElectricsAssembler 
           = mpMonodomainProblem->CreateAssembler();

        // set up initial voltage etc
        Vec voltage;        
        Vec initial_voltage = mpMonodomainProblem->CreateInitialCondition();

        // create stores of lambda, lambda_dot and old lambda
        unsigned num_quad_points = mpCardiacMechAssembler->GetTotalNumQuadPoints();

        // these are only needed if explicit
        std::vector<double> lambda;
        std::vector<double> old_lambda;
        std::vector<double> dlambda_dt;
        std::vector<NhsCellularMechanicsOdeSystem> cellmech_systems;
        EulerIvpOdeSolver euler_solver;

        // this is the active tension if explicit and the calcium conc if implicit
        std::vector<double> forcing_quantity(num_quad_points,0.0);
        
        // initial cellmechanics systems, lambda, etc, if required
        if(mUseExplicitMethod)
        {
            lambda.resize(num_quad_points, 1.0);
            old_lambda.resize(num_quad_points, 1.0);
            dlambda_dt.resize(num_quad_points, 0.0);
            cellmech_systems.resize(num_quad_points);
        }

        unsigned mech_writer_counter = 0;

        // write initial positions
        if(mWriteOutput)
        {
            OutputFileHandler output_file_handler(mOutputDirectory, true);
            out_stream p_file = output_file_handler.OpenOutputFile("results_", mech_writer_counter, ".dat");
            std::vector<Vector<double> >& deformed_position = mpCardiacMechAssembler->rGetDeformedPosition();
            for(unsigned i=0; i<deformed_position[0].size(); i++)
            {
                assert(DIM==1);
                (*p_file) << deformed_position[0](i) << "\n";
            }
        }


        unsigned counter = 0;

        TimeStepper stepper(0.0, mEndTime, mTimeStep);
        while ( !stepper.IsTimeAtEnd() )
        {
            std::cout << "**Time = " << stepper.GetTime() << "\n" << std::flush;
            
            // solve the electrics
            mpElectricsAssembler->SetTimes(stepper.GetTime(), stepper.GetNextTime(), mTimeStep);
            mpElectricsAssembler->SetInitialCondition( initial_voltage );
            voltage = mpElectricsAssembler->Solve();

            VecDestroy(initial_voltage);
            initial_voltage = voltage;
            
            // compute Ca_I at each quad point (by interpolation, using the info on which
            // electrics element the quad point is in. Then: 
            //   Explicit: Set Ca_I on the nhs systems and solve them to get the active tension
            //   Implicit: Set Ca_I on the mechanics solver
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

            // NOTE: HERE WE SHOULD REALLY CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
            // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS.. (esp for implicit)
            
            // set the active tensions
            mpCardiacMechAssembler->SetForcingQuantity(forcing_quantity);

            // solve the mechanics
            mpCardiacMechAssembler->Solve(stepper.GetTime(), stepper.GetNextTime(), stepper.GetNextTime()-stepper.GetTime());

            // if explicit store the new lambda and update lam
            if(mUseExplicitMethod)
            {
                // update lambda and dlambda_dt;
                old_lambda = lambda;
                lambda = mpCardiacMechAssembler->GetLambda();
                for(unsigned i=0; i<dlambda_dt.size(); i++)
                {
                    dlambda_dt[i] = (lambda[i] - old_lambda[i])/mTimeStep;
                }
            }

            // write
            if(mWriteOutput && (counter++)%WRITE_EVERY_NTH_TIME==0)
            {            
                OutputFileHandler output_file_handler(mOutputDirectory, false);
                out_stream p_file = output_file_handler.OpenOutputFile("results_", mech_writer_counter, ".dat");
                std::vector<Vector<double> >& deformed_position = mpCardiacMechAssembler->rGetDeformedPosition();
                for(unsigned i=0; i<deformed_position[0].size(); i++)
                {
                    assert(DIM==1);
                    (*p_file) << deformed_position[0](i) << "\n";
                }
                mech_writer_counter++;
            }                        
            
            // update the current time
            stepper.AdvanceOneTimeStep();
        }
        
        delete mpElectricsAssembler;
    }
};

#endif /*ONDEDIMCARDIACELECTROMECHANICSPROBLEM_HPP_*/
