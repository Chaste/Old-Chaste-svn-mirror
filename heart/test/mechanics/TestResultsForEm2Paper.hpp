/*

Copyright (C) University of Oxford, 2005-2009

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



#ifndef TESTRESULTSFOREM2PAPER_HPP_
#define TESTRESULTSFOREM2PAPER_HPP_


#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "ExplicitCardiacMechanicsAssembler.hpp"
#include "ImplicitCardiacMechanicsAssembler.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraturePointsGroup.hpp"
#include "NonlinearElasticityTools.hpp"
#include "ReplicatableVector.hpp"

#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


// Todo: move to pras's project folder when it lets me check that out..

class TestResultsForEm2Paper : public CxxTest::TestSuite
{
public:
    // 1. MAKE SURE  k=5 - ie change {{{return fabs(5*sin(mTime));}}} to {{{return fabs(sin(mTime));}}} in NonPhysiologicalContractionModel
    void xTestWithStretchIndependent() throw(Exception)
    {
        EXIT_IF_PARALLEL; 
        
        QuadraticMesh<2> mesh(1.0, 1.0, 5, 5);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        // NONPHYSIOL 1 - contraction model is of the form sin(t)
        ExplicitCardiacMechanicsAssembler<2> expl_solver(NONPHYSIOL1,&mesh,"TestExpVsImpStretchIndep_Exp",fixed_nodes,&law);
        ImplicitCardiacMechanicsAssembler<2> impl_solver(NONPHYSIOL1,&mesh,"TestExpVsImpStretchIndep_Imp",fixed_nodes,&law);

        expl_solver.WriteOutput(0);
        impl_solver.WriteOutput(0);

        unsigned counter = 1;

        double t0 = 0.0;
        double t1 = 10.0; 
        double dt = 0.01;

        for(double t=t0; t<t1; t+=dt)
        {
            std::cout << "\n **** t = " << t << " ****\n" << std::flush;
            
            expl_solver.SetWriteOutput(false);
            expl_solver.Solve(t,t+dt,dt);
            expl_solver.SetWriteOutput();
            expl_solver.WriteOutput(counter);

            impl_solver.SetWriteOutput(false);
            impl_solver.Solve(t,t+dt,dt);
            impl_solver.SetWriteOutput();
            impl_solver.WriteOutput(counter);

            // computations should be identical
            TS_ASSERT_EQUALS(expl_solver.GetNumNewtonIterations(), impl_solver.GetNumNewtonIterations()); 
            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](0),  impl_solver.rGetDeformedPosition()[i](0), 1e-10);
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](1),  impl_solver.rGetDeformedPosition()[i](1), 1e-10);
            }
            
            counter++;
        }
        
        /* MATLAB/OCTAVE CODE
        n = 1000;
        end_time = 10;
        nodes = 5;
        
        res = zeros(n,3);
        for i=0:n
          expl=load(['TestExpVsImpStretchIndep_Exp/solution_',num2str(i),'.nodes']);
          impl=load(['TestExpVsImpStretchIndep_Imp/solution_',num2str(i),'.nodes']);
        
          res(i+1,1) = end_time*i/n;
          res(i+1,2) = expl(nodes+1,1);
          res(i+1,3) = impl(nodes+1,1);
        end;
        
        plot(res(:,1),res(:,2),'r');
        plot(res(:,1),res(:,3),'b');
        */
    }

    
    // 0. MAKE SURE stiffness in NonPhysiologicalContractionModel is equal to 1.0
    // 1. run with dt=0.01, - will run fine
    // 2. run with dt=0.1,  
    // 3. run with dt=1.0  
    void xTestWithStretchDependentStable() throw(Exception)
    {
        double dt = 0.01;
        
        EXIT_IF_PARALLEL; 
        
        QuadraticMesh<2> mesh(1.0, 1.0, 5, 5);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        // NONPHYSIOL 2 - contraction model is of the form lam*sin(t)
        std::stringstream ss1,ss2;
        ss1 << "TestExpVsImpStretchDep_k1_dt" << dt << "_Exp";
        ss2 << "TestExpVsImpStretchDep_k1_dt" << dt << "_Imp";
        
        ExplicitCardiacMechanicsAssembler<2> expl_solver(NONPHYSIOL2,&mesh,ss1.str(),fixed_nodes,&law);
        ImplicitCardiacMechanicsAssembler<2> impl_solver(NONPHYSIOL2,&mesh,ss2.str(),fixed_nodes,&law);

        expl_solver.WriteOutput(0);
        impl_solver.WriteOutput(0);

        unsigned counter = 1;

        double t0 = 0.0;
        double t1 = 10.0; 

        for(double t=t0; t<t1; t+=dt)
        {
            std::cout << "\n **** t = " << t << " ****\n" << std::flush;
            
            expl_solver.SetWriteOutput(false);
            expl_solver.Solve(t,t+dt,dt);
            expl_solver.SetWriteOutput();
            expl_solver.WriteOutput(counter);

            impl_solver.SetWriteOutput(false);
            impl_solver.Solve(t,t+dt,dt);
            impl_solver.SetWriteOutput();
            impl_solver.WriteOutput(counter);
            
            counter++;
        }
        
        /* MATLAB/OCTAVE CODE
        dt=1;
        n = 10;
        end_time = 10;
        nodes = 5;
        
        res = zeros(n,3);
        for i=0:n
          expl=load(['TestExpVsImpStretchDep_k1_dt',num2str(dt),'_Exp/solution_',num2str(i),'.nodes']);
          impl=load(['TestExpVsImpStretchDep_k1_dt',num2str(dt),'_Imp/solution_',num2str(i),'.nodes']);
        
          res(i+1,1) = end_time*i/n;
          res(i+1,2) = expl(nodes+1,1);
          res(i+1,3) = impl(nodes+1,1);
        end;
        
        plot(res(:,1),res(:,2),'r*');
        plot(res(:,1),res(:,3),'b*');
        
        % then (outside m file) something like: save -ascii "~/Desktop/??" res
        */
    }


    // 0a. MAKE SURE stiffness in NonPhysiologicalContractionModel is equal to 5
    // 0b. make sure following line commented: //ss1 << "TestExpVsImpStretchDep_k5_h1_Exp";
    // 0c. make sure nodes=5
    // 1. run with dt=0.01
    // 2. comment out implicit solve and run with dt=0.001
    // 3. comment out implicit solve and run with dt=0.0001
    // 4. do nodes=5, dt=0.01, uncomment "//ss1 << "TestExpVsImpStretchDep_k5_h1_Exp";" and run (without implicit)
    void xTestWithStretchDependentUnstable() throw(Exception)
    {
        double dt = 0.01;
        unsigned nodes = 5;
        
        EXIT_IF_PARALLEL; 
        
        QuadraticMesh<2> mesh(1.0, 1.0, nodes, nodes);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        // NONPHYSIOL 2 - contraction model is of the form 5*lam*sin(t)
        std::stringstream ss1,ss2;

        //ss1 << "TestExpVsImpStretchDep_k5_h1_Exp";

        ss1 << "TestExpVsImpStretchDep_k5_dt" << dt << "_Exp";
        ss2 << "TestExpVsImpStretchDep_k5_dt" << dt << "_Imp";
        
        ExplicitCardiacMechanicsAssembler<2> expl_solver(NONPHYSIOL2,&mesh,ss1.str(),fixed_nodes,&law);
        ImplicitCardiacMechanicsAssembler<2> impl_solver(NONPHYSIOL2,&mesh,ss2.str(),fixed_nodes,&law);

        expl_solver.WriteOutput(0);
        impl_solver.WriteOutput(0);

        unsigned counter = 1;

        double t0 = 0.0;
        double t1 = 10.0; // 1.0 to be quite quick (min stretch ~=0.88), make this 5 say for min stretch < 0.7

        // do all implicit first as exp will crash
        for(double t=t0; t<t1; t+=dt)
        {
            std::cout << "\n **** t = " << t << " (imp) ****\n" << std::flush;

            impl_solver.SetWriteOutput(false);
            impl_solver.Solve(t,t+dt,dt);
            impl_solver.SetWriteOutput();
            impl_solver.WriteOutput(counter);
            
            counter++;
        }
            
        counter = 0;
        for(double t=t0; t<t1; t+=dt)
        {
            std::cout << "\n **** t = " << t << " (exp) ****\n" << std::flush;
            
            expl_solver.SetWriteOutput(false);
            expl_solver.Solve(t,t+dt,dt);
            expl_solver.SetWriteOutput();
            expl_solver.WriteOutput(counter);
            
            counter++;
        }
        
        /* MATLAB/OCTAVE CODE
          
         **** EXPLICT AND IMPLICIT CASE ****
        dt=0.01;
        n = 1000;
        crash_n = 119; % for k=0.01 - check before running
        
        
        end_time = 10;
        nodes = 5;
        
        res = zeros(n,3);
        for i=0:n
          impl=load(['TestExpVsImpStretchDep_k5_dt',num2str(dt),'_Imp/solution_',num2str(i),'.nodes']);
        
          res(i+1,1) = end_time*i/n;
          res(i+1,3) = impl(nodes+1,1);
        end;
        
        for i=0:crash_n
          expl=load(['TestExpVsImpStretchDep_k5_dt',num2str(dt),'_Exp/solution_',num2str(i),'.nodes']);
        
          res(i+1,2) = expl(nodes+1,1);
        end;
        
        
        plot(res(:,1),res(:,3),'b');
        plot(res(1:crash_n-1,1),res(1:crash_n-1,2),'r*');
        
        **** EXPLICT-ONLY CASE ****
        dt=0.001;
        crash_n = 805 %% check - 805 for dt=0.001, 6652 for 0.0001
        
        end_time = 10;
        n = 10/dt;
        nodes = 5;
        
        explicit_res = zeros(crash_n+1,2);        
        for i=0:crash_n
          expl=load(['TestExpVsImpStretchDep_k5_dt',num2str(dt),'_Exp/solution_',num2str(i),'.nodes']);
          explicit_res(i+1,1) = end_time*i/n;
          explicit_res(i+1,2) = expl(nodes+1,1);
        end;
        
        plot(explicit_res(:,1),explicit_res(:,2),'b');
        
        **** EXPLICIT h=1 CASE ****
        dt=0.01;
        end_time = 10;
        n = 10/dt;
        nodes = 1;
        
        explicit_res = zeros(n,2);
        for i=0:n
          expl=load(['TestExpVsImpStretchDep_k5_h1_Exp/solution_',num2str(i),'.nodes']);
          explicit_res(i+1,1) = end_time*i/n;
          explicit_res(i+1,2) = expl(nodes+1,1);
        end;
    
        plot(explicit_res(:,1),explicit_res(:,2),'b');        
        */
    } 


    // 1. edit output dir, dx and dt here, and in CardiacElectroMechanicsProblem to use
    // implicit or explicit for kerchoffs as desired
    void xTestRunWithKerchoffs() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(KERCHOFFS2003,
 						     1,   /* width (cm) */
                                                     10,  /* mech mesh size*/
                                                     100, /* elec elem each dir */
                                                     &cell_factory,
                                                     100,  /* end time */                    // dies at 7.28 with explicit
                                                     1,    /* n times 0.01ms mech dt */
                                                     0.01, /* contraction model ode timestep */
                                                     "TestExplicitWithKerchoffs_h0.1_dt001");

        c_vector<double,2> pos;
        pos(0) = 1.0;
        pos(1) = 0.0;
        
        problem.SetWatchedPosition(pos);
        problem.SetNoElectricsOutput();
        problem.Initialise();

        problem.Solve();
    }
};


#endif /*TESTRESULTSFOREM2PAPER_HPP_*/
