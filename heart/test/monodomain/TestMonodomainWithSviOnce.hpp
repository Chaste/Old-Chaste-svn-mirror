/*

Copyright (C) University of Oxford, 2005-2011

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
#include "ZeroStimulusCellFactory.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MatrixBasedMonodomainSolver.hpp"


// stimulate a block of cells (an interval in 1d, a block in a corner in 2d)
template<unsigned DIM>
class BlockCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactory()
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
        assert(DIM<3);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[0];
        double y;
        if(DIM==2)
        {
            y = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[1];
        }

        if (    (DIM==1 && fabs(x)<0.02+1e-6)
             || (DIM==2 && fabs(x)<0.02+1e-6 && fabs(y)<0.02+1e-6) )  
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

    void TestConductionVelocityInCrossFibreDirection2d() throw(Exception)
    {
 
        ReplicatableVector final_voltage_ici;
        ReplicatableVector final_voltage_svi;

        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.01);
        
        // much lower conductivity in cross-fibre direction - ICI will struggle
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.17));
        
   
        // SVI - state variable interpolation
        {
            ///\todo #1462 - TetrahedralMesh<2,2> mesh messes with Halo nodes causing GetCardiac.. to fail;
            //DistributedTetrahedralMesh<2,2> mesh;
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.04, 0.04);
            

            HeartConfig::Instance()->SetOutputDirectory("MonodomainSvi2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            BlockCellFactory<2> cell_factory;
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();

            monodomain_problem.Solve();
                
            final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());
        }
        
    }

    
};

#endif /*TESTMONODOMAINWITHSVI_HPP_*/
