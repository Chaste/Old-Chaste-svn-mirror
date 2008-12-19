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


#ifndef TESTBIDOMAINWITHBATHASSEMBLERLONG_HPP_
#define TESTBIDOMAINWITHBATHASSEMBLERLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainWithBathAssembler.hpp"
#include "TetrahedralMesh.hpp"

template<unsigned DIM>
class BathCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    SimpleStimulus* mpStimulus;
    c_vector<double,DIM> mStimulatedPoint;

public:
    BathCellFactory(double stimulusMagnitude, c_vector<double,DIM> stimulatedPoint) : AbstractCardiacCellFactory<DIM>()
    {
        // set the new stimulus
        mpStimulus = new SimpleStimulus(stimulusMagnitude, 0.5);
        mStimulatedPoint = stimulatedPoint;
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        // stimulate centre node normally.. 
        bool is_centre;
        
        if (DIM==1)
        {
            is_centre = (fabs(this->mpMesh->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6);
        }
        else if (DIM==2)
        {
            is_centre = (    (fabs(this->mpMesh->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6) 
                          && (fabs(this->mpMesh->GetNode(node)->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6) );
        }
        else
        {
            is_centre = (    (fabs(this->mpMesh->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6) 
                          && (fabs(this->mpMesh->GetNode(node)->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6) 
                          && (fabs(this->mpMesh->GetNode(node)->GetPoint()[2]-mStimulatedPoint(2)) < 1e-6) );
        }
        
        if (is_centre)
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, mpStimulus, this->mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, this->mpZeroStimulus, this->mpZeroStimulus);
        }
    }

    ~BathCellFactory(void)
    {
        delete mpStimulus;
    }
};




class TestBidomainWithBathAssemblerLong : public CxxTest::TestSuite
{
public:
    void Test3dBathIntracellularStimulation() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_3d");
                        
        c_vector<double,3> centre;
        centre(0) = 0.05;                        
        centre(1) = 0.05;
        centre(2) = 0.05;
        BathCellFactory<3> cell_factory(-2.5e7, centre); // stimulates x=0.05 node
  
        BidomainProblem<3> bidomain_problem( &cell_factory, true );

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(10,10,10);
        mesh.Scale(1.0/100.0, 1.0/100.0, 1.0/100.0);
        
        // Set everything outside a central sphere (radius 0.4) to be bath
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            double y = mesh.GetElement(i)->CalculateCentroid()[1];
            double z = mesh.GetElement(i)->CalculateCentroid()[2];
            if( sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05) + (z-0.05)*(z-0.05)) > 0.04 )
            {
                mesh.GetElement(i)->SetRegion(1);
            }
        }

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Solve();
        
        Vec sol = bidomain_problem.GetVoltage();
        ReplicatableVector sol_repl(sol);

        // test V = 0 for all bath nodes
        for(unsigned i=0; i<mesh.GetNumNodes(); i++) 
        {
            if(mesh.GetNode(i)->GetRegion()==1) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
        }
        
        // a hardcoded values
        TS_ASSERT_DELTA(sol_repl[2*404], 39.7258, 1e-3);
    }
};


#endif /*TESTBIDOMAINWITHBATHASSEMBLERLONG_HPP_*/
