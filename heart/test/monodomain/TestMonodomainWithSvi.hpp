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

#ifndef TESTMONODOMAINWITHSVI_HPP_
#define TESTMONODOMAINWITHSVI_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"


class BlockCellFactory1d : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactory1d()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        if ( fabs(x)<0.01+1e-6 )  
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

class TestMonodomainWithSvi : public CxxTest::TestSuite
{
public:
    void TestOnFineMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;
        
        double h=0.001; 
        
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(h, 1.0);

        ReplicatableVector final_voltage_nci;
        ReplicatableVector final_voltage_svi;

        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        {        
            HeartConfig::Instance()->SetOutputDirectory("MonodomainNciFine");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            BlockCellFactory1d cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
        
            monodomain_problem.Solve();
                
            final_voltage_nci.ReplicatePetscVector(monodomain_problem.GetSolution());
        }

        {        
            HeartConfig::Instance()->SetOutputDirectory("MonodomainSviFine");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            BlockCellFactory1d cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
        
            monodomain_problem.Solve();
                
            final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());
        }
        
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // visually checked they agree at this mesh resolution, and chosen tolerance from results
            TS_ASSERT_DELTA(final_voltage_nci[i], final_voltage_svi[i], 0.3);
            
            if(final_voltage_nci[i]>-80)
            {
                // shouldn't be exactly equal, as long as away from resting potential
                TS_ASSERT_DIFFERS(final_voltage_nci[i], final_voltage_svi[i]);
            }
        }
    }
};

#endif /*TESTMONODOMAINWITHSVI_HPP_*/
