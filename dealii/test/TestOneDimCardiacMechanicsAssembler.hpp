#ifndef TEST1DCARDIACMECHANICSASSEMBLER_HPP_
#define TEST1DCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "OneDimCardiacMechanicsAssembler.hpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "FiniteElasticityTools.hpp"
#include "DealiiMonodomainAssembler.hpp" // to be replaced by [Mono/Bi]Dg0Assmebler...
#include "InitialStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AbstractCardiacCell.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "NHSCellularMechanicsOdeSystem.hpp"

class TestOneDimCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestPoleZero3dIn1dLaw()
    {
        PoleZero3dIn1dLaw law;

        TS_ASSERT_DELTA( law.GetT(0), 0.0, 1e-12 );
        TS_ASSERT_DELTA( law.GetT(0.1), 2.0809, 1e-3 );
        TS_ASSERT_DELTA( law.GetT(0.2), 8.5158, 1e-3 );
        TS_ASSERT_DELTA( law.GetT(0.3), 37.0291, 1e-3 );

        TS_ASSERT_DELTA( law.GetT(-0.1), -0.5023, 1e-3 );
        TS_ASSERT_DELTA( law.GetT(-0.2), -2.2589, 1e-3 );

        law.SetUpStores();

        TS_ASSERT_DELTA( law.GetT(-0.1), -0.5023, 1e-3 );
        TS_ASSERT_DELTA( law.GetT(-0.2), -2.2589, 1e-3 );
    }


    void TestSimple()
    {
        Triangulation<1> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(7);
        
        OneDimCardiacMechanicsAssembler mechanics(&mesh);

        std::vector<double> active_tension(mechanics.GetTotalNumQuadPoints(), 1);
        
        mechanics.SetActiveTension( active_tension );
        mechanics.Solve();
        
        std::vector<Vector<double> > undeformed_position = mechanics.rGetUndeformedPosition();
        std::vector<Vector<double> > deformed_position = mechanics.rGetDeformedPosition();
        
        // not i=0 as X=0
        double factor = deformed_position[0](1)/undeformed_position[0](1);
        
        TS_ASSERT_LESS_THAN(factor, 1.0);
        
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            if(undeformed_position[0](i)!=0)
            {
                TS_ASSERT_DELTA(deformed_position[0](i)/undeformed_position[0](i), factor, 1e-4);
            }
        }
    }


    void TestOneDimElectroMechanics()
    {
        // create a single mesh for both electrics and mechanics
        Triangulation<1> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 10);
        mesh.refine_global(7);
        
        // create two stimulus functions
        AbstractStimulusFunction* p_stimulus = new InitialStimulus(-1e6, 0.5);;
        AbstractStimulusFunction* p_zero_stim = new ZeroStimulus;
        EulerIvpOdeSolver euler_solver;
        
        double time_step = 0.1;
        
        // create the cells, stimulating one face
        std::vector< AbstractCardiacCell* > cells(mesh.n_vertices());
        
        TriangulationVertexIterator<1> vertex_iter(&mesh);
        while(!vertex_iter.ReachedEnd())
        {
            if(vertex_iter.GetVertex()[0]==0)
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new BackwardEulerLuoRudyIModel1991(time_step, p_stimulus);
            }
            else
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new BackwardEulerLuoRudyIModel1991(time_step, p_zero_stim);
            }
            
            vertex_iter.Next();
        }
        unsigned Ca_i_index = cells[0]->GetStateVariableNumberByName("CaI");
                
        // create the dealii monodomain assembler
        DealiiMonodomainAssembler<1> monodomain_assembler(&mesh, cells);


        // create an assembler for the mechanics
        OneDimCardiacMechanicsAssembler cardiac_mech_assembler(&mesh);


        // create stores of lambda, lambda_dot and old lambda
        unsigned num_quad_points = cardiac_mech_assembler.GetTotalNumQuadPoints();
        std::vector<double> lambda(num_quad_points, 1.0);
        std::vector<double> old_lambda(num_quad_points, 1.0);
        std::vector<double> dlambda_dt(num_quad_points, 0.0);
        
        std::vector<double> active_tension(num_quad_points, 0.0);

        // create NHS systems for each quad point in the mesh
        std::vector<NHSCellularMechanicsOdeSystem> cellmech_systems(num_quad_points);
        
        double time = 0;
        double end_time = 1000; 
        double dt = time_step;
        
      //  double time_since_last_mech_solve = 0.0;
        unsigned mech_writer_counter = 0;
     //   double mech_solving_time = 0.1;    // solve the mechanics every 0.1 ms
        

        OutputFileHandler output_file_handler("OneDimElectroMech", true);
        std::stringstream ss;
        ss << "results_" << mech_writer_counter << ".dat"; 
        out_stream p_file = output_file_handler.OpenOutputFile(ss.str());
        std::vector<Vector<double> >& deformed_position = cardiac_mech_assembler.rGetDeformedPosition();
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            (*p_file) << deformed_position[0](i) << "\n";
        }

        while(time < end_time)
        {
            std::cout << "Time = " << time << "\n";
     
            // solve for the voltage (with solves the electrical cells)
            monodomain_assembler.Solve(time, time+dt, 1);
                
//            std::cout << "Voltage:\n";
//            for(uint i=0; i<cells.size(); i++)
//            {
//                std::cout << monodomain_assembler.rGetVoltage()[i] << " ";
//            }
//            std::cout << "\n";
//
//            std::cout << "Ca_I\n";
//            for(uint i=0; i<cells.size(); i++)
//            {
//                std::cout << cells[i]->rGetStateVariables()[Ca_i_index] << " "; 
//            }
//            std::cout << "\n";

//                std::cout << "average_Ca_I:\n";
            // set Ca_I on each NHS system. Set lamda and dlamda_dt at the same time.                

            unsigned current_quad_index = 0;                
            for(Triangulation<1>::cell_iterator element_iter = mesh.begin_active(); 
                element_iter!=mesh.end();
                element_iter++)
            {
                // average the Ca_I values from the nodes of this element
                double average_Ca_I = 0;
                for(unsigned i=0; i<GeometryInfo<1>::vertices_per_cell; i++)
                {
                    unsigned node_index = element_iter->vertex_index(i);
                    average_Ca_I += cells[node_index]->rGetStateVariables()[Ca_i_index];
                }
                average_Ca_I /= GeometryInfo<1>::vertices_per_cell;
                
//                std::cout << average_Ca_I << " ";
                
                // apply to all nhs systems for quad points in this element
                for(unsigned i=0; i<cardiac_mech_assembler.GetNumQuadPointsPerElement(); i++)
                {
                    cellmech_systems[current_quad_index].SetLambda1DerivativeAndCalciumI(lambda[current_quad_index], 
                                                                                         dlambda_dt[current_quad_index],
                                                                                         average_Ca_I);
                    current_quad_index++;                                                                                         
                }
            }

//            std::cout << "\nActive Tension:\n";
            // solve the cellular mechanics model and get the active tension
            for(unsigned i=0; i<cellmech_systems.size(); i++)
            {
                euler_solver.SolveAndUpdateStateVariable(&cellmech_systems[i], time, time+dt, dt);
                active_tension[i] = cellmech_systems[i].GetActiveTension();
            }
            
            // NOTE: HERE WE SHOULD REALLY CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
            // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS..
            
            // set the active tensions
            cardiac_mech_assembler.SetActiveTension(active_tension);
            
            
//            if(time_since_last_mech_solve > mech_solving_time)
  //          {
                // solve the mechanics system
                cardiac_mech_assembler.Solve();
                
                OutputFileHandler output_file_handler("OneDimElectroMech", false);
                std::stringstream ss;
                ss << "results_" << mech_writer_counter << ".dat"; 
                out_stream p_file = output_file_handler.OpenOutputFile(ss.str());
                for(unsigned i=0; i<deformed_position[0].size(); i++)
                {
                    (*p_file) << deformed_position[0](i) << "\n";
                }
                
                mech_writer_counter++;
               // time_since_last_mech_solve = 0.0;
            //}

            // update lambda and dlambda_dt;
            old_lambda = lambda;
            lambda = cardiac_mech_assembler.GetLambda();
            for(unsigned i=0; i<dlambda_dt.size(); i++)
            {
                dlambda_dt[i] = (lambda[i] - old_lambda[i])/dt;
            }
            
            std::vector<Vector<double> >& deformed_position = cardiac_mech_assembler.rGetDeformedPosition();
            
    //        time_since_last_mech_solve += dt;
            time += dt;
        }
        
//        // write the file solution, if it hasn't just been written.        
//        if(time_since_last_mech_solve > 0.0)
//        {   
//        }
        
    }
};
#endif /*TEST1DCARDIACMECHANICSASSEMBLER_HPP_*/
