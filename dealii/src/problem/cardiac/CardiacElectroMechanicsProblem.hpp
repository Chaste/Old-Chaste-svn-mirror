#ifndef CARDIACELECTROMECHANICSPROBLEM_HPP_
#define CARDIACELECTROMECHANICSPROBLEM_HPP_

#include "MonodomainProblem.hpp"
#include "OneDimCardiacMechanicsAssembler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "FiniteElasticityTools.hpp"
#include "NHSCellularMechanicsOdeSystem.hpp"


/* todos:
 * 
 * add comments
 * add tests
 * 
 * tidy
 * make two-mesh
 * move mesh stuff out
 * don't output every timestep
 */


///\todo -make this struct
template<unsigned DIM>
class ElementAndWeights
{
public:
    unsigned ElementNum;
    c_vector<double,DIM+1> Weights;  
};


template<unsigned DIM>
class CardiacElectroMechanicsProblem
{
private: 
    MonodomainProblem<DIM>* mpMonodomainProblem;
    OneDimCardiacMechanicsAssembler* mpCardiacMechAssembler;

    double mEndTime;
    double mTimeStep;    
    
    ConformingTetrahedralMesh<DIM,DIM>* mpElectricsMesh;
    Triangulation<DIM>*                 mpMechanicsMesh;

    std::vector<ElementAndWeights<DIM> > mElementAndWeightsForQuadPoints;

    std::string mOutputDirectory;
        
public:
    CardiacElectroMechanicsProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   double endTime,
                                   double timeStep,
                                   std::string outputDirectory)
    {
        assert(pCellFactory != NULL);
        mpMonodomainProblem = new MonodomainProblem<DIM>(pCellFactory);
        
        assert(endTime > 0);
        mEndTime = endTime;
        mTimeStep = timeStep;
        
        assert(outputDirectory!="");
        mOutputDirectory = outputDirectory;
                
        mpElectricsMesh = new ConformingTetrahedralMesh<DIM,DIM>();
        unsigned num_elem = 128;
        mpElectricsMesh->ConstructLinearMesh(num_elem);
        mpElectricsMesh->Scale(1.0/num_elem);
        
        assert(DIM==1);
        mpMonodomainProblem->SetMesh(mpElectricsMesh);
        mpMonodomainProblem->Initialise();
        
        mpMechanicsMesh = new Triangulation<DIM>();
        GridGenerator::hyper_cube(*mpMechanicsMesh, 0.0, 1);
        mpMechanicsMesh->refine_global(7);
        
        assert(mpMechanicsMesh->n_vertices()==mpElectricsMesh->GetNumNodes());

        mpCardiacMechAssembler = new OneDimCardiacMechanicsAssembler(mpMechanicsMesh);
        
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
    
    
    void Solve()
    {
        AbstractDynamicAssemblerMixin<DIM,DIM,1>* mpElectricsAssembler 
           = mpMonodomainProblem->CreateAssembler();
        Vec voltage;        
        Vec initial_voltage = mpMonodomainProblem->CreateInitialCondition();

        // create stores of lambda, lambda_dot and old lambda
        unsigned num_quad_points = mpCardiacMechAssembler->GetTotalNumQuadPoints();
        std::vector<double> lambda(num_quad_points, 1.0);
        std::vector<double> old_lambda(num_quad_points, 1.0);
        std::vector<double> dlambda_dt(num_quad_points, 0.0);
        std::vector<double> active_tension(num_quad_points, 0.0);

        // create NHS systems for each quad point in the mesh
        std::vector<NHSCellularMechanicsOdeSystem> cellmech_systems(num_quad_points);

        EulerIvpOdeSolver euler_solver;

        unsigned mech_writer_counter = 0;

        OutputFileHandler output_file_handler(mOutputDirectory, true);
        out_stream p_file = output_file_handler.OpenOutputFile("results_", mech_writer_counter, ".dat");
        std::vector<Vector<double> >& deformed_position = mpCardiacMechAssembler->rGetDeformedPosition();
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            assert(DIM==1);
            (*p_file) << deformed_position[0](i) << "\n";
        }


unsigned counter = 0;

        TimeStepper stepper(0.0, mEndTime, mTimeStep);
        while ( !stepper.IsTimeAtEnd() )
        {
            mpElectricsAssembler->SetTimes(stepper.GetTime(), stepper.GetNextTime(), mTimeStep);
            mpElectricsAssembler->SetInitialCondition( initial_voltage );

            voltage = mpElectricsAssembler->Solve();

            VecDestroy(initial_voltage);
            initial_voltage = voltage;
            
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

                cellmech_systems[i].SetLambda1DerivativeAndCalciumI(lambda[i], 
                                                                    dlambda_dt[i],
                                                                    interpolated_Ca_I);
            }

            // solve the cellular mechanics model and get the active tension
            for(unsigned i=0; i<cellmech_systems.size(); i++)
            {
                euler_solver.SolveAndUpdateStateVariable(&cellmech_systems[i], stepper.GetTime(), stepper.GetNextTime(), mTimeStep);
                active_tension[i] = cellmech_systems[i].GetActiveTension();
            }
            
            // NOTE: HERE WE SHOULD REALLY CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
            // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS..
            
            // set the active tensions
            mpCardiacMechAssembler->SetActiveTension(active_tension);

            mpCardiacMechAssembler->Solve();

            // update lambda and dlambda_dt;
            old_lambda = lambda;
            lambda = mpCardiacMechAssembler->GetLambda();
            for(unsigned i=0; i<dlambda_dt.size(); i++)
            {
                dlambda_dt[i] = (lambda[i] - old_lambda[i])/mTimeStep;
            }

if((counter++)%1==0)
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

#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
