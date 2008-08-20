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


#ifndef TESTHEARTCONFIG_HPP_
#define TESTHEARTCONFIG_HPP_

#include <cxxtest/TestSuite.h>
#include "HeartConfig.hpp"

class TestHeartConfig : public CxxTest::TestSuite
{
public :
    void TestHeartConfigBasic()
    {
        double chi = HeartConfig::Instance()->DefaultParameters()->Physiological().SurfaceAreaToVolumeRatio().get();
        TS_ASSERT_EQUALS(chi, 1400);

        HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteParametersFullFormat.xml");

        chi = HeartConfig::Instance()->UserParameters()->Physiological().SurfaceAreaToVolumeRatio().get();
        TS_ASSERT_EQUALS(chi, 1400);

        double capacitance = HeartConfig::Instance()->UserParameters()->Physiological().Capacitance().get();
        TS_ASSERT_EQUALS(capacitance, 1.0);

        double conductivity_1 = HeartConfig::Instance()->UserParameters()->Physiological().IntracellularConductivities().get().longi();
        double conductivity_2 = HeartConfig::Instance()->UserParameters()->Physiological().ExtracellularConductivities().get().longi();

        TS_ASSERT_EQUALS(conductivity_1, 1.75);
        TS_ASSERT_EQUALS(conductivity_2, 7.0);
    }

    void TestUserProvidedDifferentFromDefault()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteParametersFullFormat.xml");

        ionic_model_type default_ionic_model = HeartConfig::Instance()->DefaultParameters()->Simulation().IonicModel().get();
        TS_ASSERT_EQUALS(default_ionic_model, ionic_model_type::LuoRudyIModel1991OdeSystem);

        ionic_model_type user_ionic_model = HeartConfig::Instance()->UserParameters()->Simulation().IonicModel().get();
        TS_ASSERT_EQUALS(user_ionic_model, ionic_model_type::FaberRudy2000Version3);

        ionic_model_type get_ionic_model = HeartConfig::Instance()->GetIonicModel();
        TS_ASSERT_EQUALS(user_ionic_model, get_ionic_model);

     }

    void TestGetFunctions()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteParametersFullFormat.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(),
                         10.0);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDomain(),
                         domain_type::Mono);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetIonicModel(),
                          ionic_model_type::FaberRudy2000Version3);

        TS_ASSERT(HeartConfig::Instance()->GetCreateSlab());
        TS_ASSERT(!HeartConfig::Instance()->GetLoadMesh());

        c_vector<double, 3> slab_dimensions;
        HeartConfig::Instance()->GetSlabDimensions(slab_dimensions);
        
        TS_ASSERT_EQUALS(slab_dimensions[0], 4.0);
        TS_ASSERT_EQUALS(slab_dimensions[1], 0.1);
        TS_ASSERT_EQUALS(slab_dimensions[2], 2.0);
        
        double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();
        TS_ASSERT_EQUALS(inter_node_space, 0.1);

        std::vector<SimpleStimulus> stimuli_applied;
        std::vector<ChasteCuboid> stimulated_areas;
        HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_areas);

        TS_ASSERT_EQUALS(stimuli_applied.size(), 2u);
        TS_ASSERT_EQUALS(stimulated_areas.size(), 2u);

        TS_ASSERT_EQUALS(stimuli_applied[0].GetStimulus(0), -25500.0);
        TS_ASSERT_EQUALS(stimuli_applied[0].GetStimulus(0.6), 0.0);

        TS_ASSERT(stimulated_areas[1].DoesContain(ChastePoint<3>(-2, 0, -2)));
        TS_ASSERT( ! stimulated_areas[1].DoesContain(ChastePoint<3>(-6, -6, -6)));

        std::vector<ChasteCuboid> cell_heterogeneity_areas;
        std::vector<double> scale_factor_gks;
        std::vector<double> scale_factor_ito;
        HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                        scale_factor_gks,
                                                        scale_factor_ito);

        TS_ASSERT(cell_heterogeneity_areas[0].DoesContain(ChastePoint<3>(-1.0, 0, 0)));
        TS_ASSERT_EQUALS(scale_factor_gks[1], 1.154);
        TS_ASSERT_EQUALS(scale_factor_ito[2], 1);

        std::vector<ChasteCuboid> conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        TS_ASSERT(conductivities_heterogeneity_areas[0].DoesContain(ChastePoint<3>(1.95, 0, 0)));
        TS_ASSERT_EQUALS(intra_h_conductivities[0][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities[0][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities[1][0], 0.75);
        
        c_vector<double, 3> extra_conductivities;
        HeartConfig::Instance()->GetExtracellularConductivities(extra_conductivities);
        TS_ASSERT_EQUALS(extra_h_conductivities[1][0], extra_conductivities[0]);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "ChasteResults");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputFilenamePrefix(), "SimulationResults");

        c_vector<double, 3> intra_conductivities;
        HeartConfig::Instance()->GetIntracellularConductivities(intra_conductivities);
        TS_ASSERT_EQUALS(intra_conductivities[0], 1.75);
        TS_ASSERT_EQUALS(intra_conductivities[1], 1.75);
        TS_ASSERT_EQUALS(intra_conductivities[2], 1.75);

        TS_ASSERT_EQUALS(extra_conductivities[0], 7.0);
        TS_ASSERT_EQUALS(extra_conductivities[1], 7.0);
        TS_ASSERT_EQUALS(extra_conductivities[2], 7.0);
        
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio(), 1400.0);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCapacitance(), 1.0);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 0.025);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 0.025);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetPdeTimeStep(), 0.05);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetPrintingTimeStep(), 1.0);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseAbsoluteTolerance(), false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseRelativeTolerance(), true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetRelativeTolerance(), 1e-6);

        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "gmres")==0);
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "ilu")==0);

        HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteParametersLoadMesh.xml");

        TS_ASSERT(!HeartConfig::Instance()->GetCreateSlab());
        TS_ASSERT(HeartConfig::Instance()->GetLoadMesh());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshName(), "foo");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetConductivityMedia(), media_type::NoFibreOrientation);

        //Try reading through an empty parameters file into a fully-populated default
        HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteEmpty.xml");
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/ChasteParametersFullFormat.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        HeartConfig::Instance()->Reset();
    }    
    
    void TestGetIsMeshProvided()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteEmpty.xml");
		TS_ASSERT_EQUALS(HeartConfig::Instance()->GetIsMeshProvided(), false);
		TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->GetLoadMesh())
		TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->GetCreateSlab())
		
		HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteParametersLoadMesh.xml");
		TS_ASSERT_EQUALS(HeartConfig::Instance()->GetIsMeshProvided(), true);        
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetLoadMesh(), true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCreateSlab(), false);    	
    }

    void TestSetFunctions() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSimulationDuration(35.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 35.0);

        HeartConfig::Instance()->SetDomain(domain_type::Bi);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDomain(), domain_type::Bi);

        HeartConfig::Instance()->SetIonicModel(ionic_model_type::BackwardEulerLuoRudyIModel1991);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetIonicModel(), ionic_model_type::BackwardEulerLuoRudyIModel1991);

        HeartConfig::Instance()->SetOutputDirectory("NewOuputDirectory");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "NewOuputDirectory");

        HeartConfig::Instance()->SetOutputFilenamePrefix("NewSimulation");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputFilenamePrefix(), "NewSimulation");

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(-6.0, -5.0, -4.0));
        
        c_vector<double, 3> intra;
        HeartConfig::Instance()->GetIntracellularConductivities(intra);
        TS_ASSERT_EQUALS(intra[0], -6.0);
        TS_ASSERT_EQUALS(intra[1], -5.0);
        TS_ASSERT_EQUALS(intra[2], -4.0);

        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(-3.0, -2.0, -1.0));
        c_vector<double, 3> extra;
        HeartConfig::Instance()->GetExtracellularConductivities(extra);
        TS_ASSERT_EQUALS(extra[0], -3.0);
        TS_ASSERT_EQUALS(extra[1], -2.0);
        TS_ASSERT_EQUALS(extra[2], -1.0);

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio(), 2000);

        HeartConfig::Instance()->SetCapacitance(2.3);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCapacitance(), 2.3);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(1.1,2.2,4.4);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 1.1);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetPdeTimeStep(), 2.2);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetPrintingTimeStep(), 4.4);

        HeartConfig::Instance()->SetOdeTimeStep(0.1);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 0.1);

        HeartConfig::Instance()->SetPdeTimeStep(0.2);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetPdeTimeStep(), 0.2);

        HeartConfig::Instance()->SetPrintingTimeStep(0.4);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetPrintingTimeStep(), 0.4);
        
        // Test code to check consistency among TimeSteps
        // throws because argument is negative
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,0.1);
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetOdeTimeStep(0.2));
        
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(-0.1, 0.1, 0.1));
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, -0.1, 0.1));
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, -0.1));

        //Throws when we try to print more often than the pde time step
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.1));
         //Throws when printing step is not a multiple of pde time step
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.3));
        
		// Throws because ode time step is bigger than pde time step
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.3));
        

        HeartConfig::Instance()->SetUseRelativeTolerance(1e-4);
        TS_ASSERT(HeartConfig::Instance()->GetUseRelativeTolerance());
        TS_ASSERT(HeartConfig::Instance()->GetUseAbsoluteTolerance() == false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetRelativeTolerance(), 1e-4);

        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-11);
        TS_ASSERT(HeartConfig::Instance()->GetUseAbsoluteTolerance());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetAbsoluteTolerance(), 1e-11);
        //Check that relative tolerance is disabled
        TS_ASSERT(HeartConfig::Instance()->GetUseRelativeTolerance() == false);
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->GetRelativeTolerance());
        
       

        HeartConfig::Instance()->SetUseRelativeTolerance(1e-4);
        TS_ASSERT(HeartConfig::Instance()->GetUseRelativeTolerance());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetRelativeTolerance(), 1e-4);
        //Check that absolute tolerance is disabled
        TS_ASSERT(HeartConfig::Instance()->GetUseAbsoluteTolerance() == false);
        TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->GetAbsoluteTolerance());
 
        HeartConfig::Instance()->SetKSPSolver("cg");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "cg")==0);

        HeartConfig::Instance()->SetKSPSolver("gmres");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "gmres")==0);

        HeartConfig::Instance()->SetKSPSolver("symmlq");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "symmlq")==0);

		TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetKSPSolver("foobar"));

        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "jacobi")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "bjacobi")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("ilu");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "ilu")==0);
        
		TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetKSPPreconditioner("foobar"));
        

     }
    
   void TestExceptions() throw (Exception)
   { 
       //We might want to remove SetDefaultsFile()
       HeartConfig::Instance()->Reset();
       TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetDefaultsFile("DoesNotExist.xml")); 
       TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetDefaultsFile("heart/test/data/ChasteInconsistent.xml")); 
       //TS_ASSERT_THROWS_ANYTHING(HeartConfig::Instance()->SetParametersFile("heart/test/data/ChasteInconsistent.xml")); 
    }
};

#endif /*TESTHEARTCONFIG_HPP_*/
