#ifndef TESTCARDIACELECTROMECHANICS_HPP_
#define TESTCARDIACELECTROMECHANICS_HPP_


#include <cxxtest/TestSuite.h>
#include "CardiacMechanicsAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "FiniteElasticityTools.hpp"
#include "DealiiMonodomainAssembler.hpp" // to be replaced by [Mono/Bi]Dg0Assmebler...
#include "InitialStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AbstractCardiacCell.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "CardiacMechanicsAssembler.cpp"
#include "NHSCellularMechanicsOdeSystem.hpp"

class TestCardiacElectroMechanics : public CxxTest::TestSuite
{
public:
    // basic test of the dealii monodomain assembler (which is only going to be used
    // until dealii can be connected to chaste)
    void TestMonodomainDealiiAssembler() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 0.1);
        mesh.refine_global(4);
        
        AbstractStimulusFunction* p_stimulus = new InitialStimulus(-1e6, 0.5);;
        AbstractStimulusFunction* p_zero_stim = new ZeroStimulus;
        EulerIvpOdeSolver euler_solver;
        
        double time_step = 0.01;
        
        std::vector< AbstractCardiacCell* > cells(mesh.n_vertices());
        
        TriangulationVertexIterator<2> vertex_iter(&mesh);
        while(!vertex_iter.ReachedEnd())
        {
            if(vertex_iter.GetVertex()[0]==0)
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new LuoRudyIModel1991OdeSystem(&euler_solver, time_step, p_stimulus);
            }
            else
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new LuoRudyIModel1991OdeSystem(&euler_solver, time_step, p_zero_stim);
            }
            
            vertex_iter.Next();
        }
        
        DealiiMonodomainAssembler<2> monodomain_assembler(&mesh, cells);
        
        double end_time = 2;
        unsigned num_timesteps = (unsigned)floor(end_time/time_step+0.5);
        monodomain_assembler.Solve(0, end_time, num_timesteps);
        
        for(unsigned i=0; i<cells.size(); i++)
        {
            delete cells[i];
        }
        
        TS_ASSERT_DELTA(monodomain_assembler.GetVoltage()[1], 41.8787, 1e-3);
    }
    

    /************************************************************************************
     * 
     *  ElectroMechanics:
     *  
     *  Test putting together DealiiMonodomain, NHS systems and CardiacMechanics....
     * 
     *  Solves on same mesh, and Ca_I is not properly interpolated
     * 
     *  todo: Choose correct parameters, esp material law
     *        Check if supplying material law is working correctly
     *        check starting guesses for newtons/gmres
     *        Turn test into source and test properly..
     * 
     *  big todos: use two meshes and proper interpolation
     *             transfer Ca_Trop back from coarse quad points to fine mesh?
     *             alter voltage diffusion equation to take into accout deformation
     * 
     ************************************************************************************/
     
    void TestCardiacElectroMechanicsOneMesh()
    {
        // create a single mesh for both electrics and mechanics
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 0.1);
        mesh.refine_global(3);
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        // create two stimulus functions
        AbstractStimulusFunction* p_stimulus = new InitialStimulus(-1e6, 0.5);;
        AbstractStimulusFunction* p_zero_stim = new ZeroStimulus;
        EulerIvpOdeSolver euler_solver;
        
        double time_step = 0.01;
        
        // create the cells, stimulating one face
        std::vector< AbstractCardiacCell* > cells(mesh.n_vertices());
        
        TriangulationVertexIterator<2> vertex_iter(&mesh);
        while(!vertex_iter.ReachedEnd())
        {
            if(vertex_iter.GetVertex()[0]==0)
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new LuoRudyIModel1991OdeSystem(&euler_solver, time_step, p_stimulus);
            }
            else
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new LuoRudyIModel1991OdeSystem(&euler_solver, time_step, p_zero_stim);
            }
            
            vertex_iter.Next();
        }
        unsigned Ca_i_index = cells[0]->GetStateVariableNumberByName("CaI");
                
        // create the dealii monodomain assembler
        DealiiMonodomainAssembler<2> monodomain_assembler(&mesh, cells);

        // create a material law, scaling the material params to help gmres convergence
        double mat_law_scaling = 1000;

// what should this be?
        MooneyRivlinMaterialLaw<2> material_law(20/mat_law_scaling);

        // create an assembler for the mechanics
        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, 
                                                            "CardiacElectroMech/OneMesh",
                                                            &material_law);

        // create stores of lambda, lambda_dot and old lambda
        unsigned num_quad_points = cardiac_mech_assembler.GetTotalNumQuadPoints();
        std::vector<double> lambda(num_quad_points, 1.0);
        std::vector<double> old_lambda(num_quad_points, 1.0);
        std::vector<double> dlambda_dt(num_quad_points, 0.0);
        
        // the active tension, scaled by some param (see below). has to be scaled as 
        // the material law will be a scaled material law to help gmres converge
        std::vector<double> scaled_active_tension(num_quad_points, 0.0);

        // create NHS systems for each quad point in the mesh
        std::vector<NHSCellularMechanicsOdeSystem> cellmech_systems(num_quad_points);
        
        double time = 0;
        double end_time = 1;
        double dt = 0.01;
        while(time < end_time)
        {
            //std::cout << time << "\n";
     
            // solve for the voltage (with solves the electrical cells)
            monodomain_assembler.Solve(time, time+dt, 1);
                
            // set Ca_I on each NHS system. Set lamda and dlamda_dt at the same time.                
            unsigned current_quad_index = 0;                
            for(Triangulation<2>::cell_iterator element_iter = mesh.begin_active(); 
                element_iter!=mesh.end();
                element_iter++)
            {
                // average the Ca_I values from the nodes of this element
                double average_Ca_I = 0;
                for(unsigned i=0; i<GeometryInfo<2>::vertices_per_cell; i++)
                {
                    unsigned node_index = element_iter->vertex_index(i);
                    average_Ca_I += cells[node_index]->rGetStateVariables()[Ca_i_index];
                }
                average_Ca_I /= GeometryInfo<2>::vertices_per_cell;
                
                // apply to all nhs systems for quad points in this element
                for(unsigned i=0; i<cardiac_mech_assembler.GetNumQuadPointsPerElement(); i++)
                {
                    cellmech_systems[current_quad_index].SetLambda1DerivativeAndCalciumI(lambda[current_quad_index], 
                                                                                         dlambda_dt[current_quad_index],
                                                                                         average_Ca_I);
                    current_quad_index++;                                                                                         
                }
            }

            // solve the cellular mechanics model and get the active tension
            for(unsigned i=0; i<cellmech_systems.size(); i++)
            {
                euler_solver.SolveAndUpdateStateVariable(&cellmech_systems[i], time, time+dt, dt);
                
                //std::cout << cellmech_systems[i].GetActiveTension() << " ";
                
                scaled_active_tension[i] = cellmech_systems[i].GetActiveTension() / mat_law_scaling;
            }
            
            // NOTE: HERE WE SHOULD REALLY CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
            // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS..
            
            // set the active tensions
            cardiac_mech_assembler.SetActiveTension(scaled_active_tension);
            
            // solve the mechanics system
            cardiac_mech_assembler.Solve();

            // update lambda and dlambda_dt;
            old_lambda = lambda;
            lambda = cardiac_mech_assembler.GetLambda();
            for(unsigned i=0; i<dlambda_dt.size(); i++)
            {
                dlambda_dt[i] = (lambda[i] - old_lambda[i])/dt;
            }

            time += dt;
        }
    }
};
#endif /*TESTCARDIACELECTROMECHANICS_HPP_*/
