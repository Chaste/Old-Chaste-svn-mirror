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


#ifndef TESTHEARTCONFIG_HPP_
#define TESTHEARTCONFIG_HPP_

#include "UblasCustomFunctions.hpp"

#include <sys/stat.h>
#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "HeartConfig.hpp"
#include "OutputFileHandler.hpp"
#include "ChasteCuboid.hpp"
#include "Version.hpp"
#include "TetrahedralMesh.hpp"

class TestHeartConfig : public CxxTest::TestSuite
{
private:
    void setUp()
    {
        HeartConfig::Reset();
    }
public :
    void TestHeartConfigBasic()
    {
        double chi = HeartConfig::Instance()->mpDefaultParameters->Physiological().SurfaceAreaToVolumeRatio().get();
        TS_ASSERT_EQUALS(chi, 1400);

        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");

        chi = HeartConfig::Instance()->mpUserParameters->Physiological().SurfaceAreaToVolumeRatio().get();
        TS_ASSERT_EQUALS(chi, 1400);

        double capacitance = HeartConfig::Instance()->mpUserParameters->Physiological().Capacitance().get();
        TS_ASSERT_EQUALS(capacitance, 1.0);

        double conductivity_1 = HeartConfig::Instance()->mpUserParameters->Physiological().IntracellularConductivities().get().longi();
        double conductivity_2 = HeartConfig::Instance()->mpUserParameters->Physiological().ExtracellularConductivities().get().longi();

        TS_ASSERT_EQUALS(conductivity_1, 1.75);
        TS_ASSERT_EQUALS(conductivity_2, 7.0);
    }

    void TestUserProvidedDifferentFromDefault()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");

        TS_ASSERT(HeartConfig::Instance()->mpDefaultParameters->Simulation().present());
        cp::simulation_type default_sim_elt = HeartConfig::Instance()->mpDefaultParameters->Simulation().get();

        TS_ASSERT(default_sim_elt.IonicModels().present());
        TS_ASSERT(default_sim_elt.IonicModels().get().Default().Hardcoded().present());
        cp::ionic_models_available_type default_ionic_model = default_sim_elt.IonicModels().get().Default().Hardcoded().get();
        TS_ASSERT_EQUALS(default_ionic_model, cp::ionic_models_available_type::LuoRudyI);

        TS_ASSERT(HeartConfig::Instance()->mpUserParameters->Simulation().present());
        cp::simulation_type user_sim_elt = HeartConfig::Instance()->mpUserParameters->Simulation().get();
        TS_ASSERT(user_sim_elt.IonicModels().present());
        TS_ASSERT(user_sim_elt.IonicModels().get().Default().Hardcoded().present());
        cp::ionic_models_available_type user_ionic_model = user_sim_elt.IonicModels().get().Default().Hardcoded().get();
        TS_ASSERT_EQUALS(user_ionic_model, cp::ionic_models_available_type::FaberRudy2000);

        cp::ionic_models_available_type get_ionic_model = HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get();
        TS_ASSERT_EQUALS(user_ionic_model, get_ionic_model);

     }

    void TestGetFunctions()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSpaceDimension(),
                         3u);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(),
                         10.0);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDomain(),
                         cp::domain_type::Mono);

        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(),
                         cp::ionic_models_available_type::FaberRudy2000);
        
        
        TS_ASSERT(HeartConfig::Instance()->GetCreateMesh());
        TS_ASSERT(!HeartConfig::Instance()->GetLoadMesh());

        c_vector<double, 3> slab_dimensions;
        HeartConfig::Instance()->GetSlabDimensions(slab_dimensions);

        TS_ASSERT_EQUALS(slab_dimensions[0], 4.0);
        TS_ASSERT_EQUALS(slab_dimensions[1], 0.1);
        TS_ASSERT_EQUALS(slab_dimensions[2], 2.0);

        double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();
        TS_ASSERT_EQUALS(inter_node_space, 0.1);

        std::vector<boost::shared_ptr<SimpleStimulus> > stimuli_applied;
        std::vector<ChasteCuboid<3> > stimulated_areas;
        HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_areas);

        TS_ASSERT_EQUALS(stimuli_applied.size(), 2u);
        TS_ASSERT_EQUALS(stimulated_areas.size(), 2u);

        TS_ASSERT_EQUALS(stimuli_applied[0]->GetStimulus(0), -25500.0);
        TS_ASSERT_EQUALS(stimuli_applied[0]->GetStimulus(0.6), 0.0);
        
        //covering the 2D case
        std::vector<ChasteCuboid<2> > stimulated_areas_2D;
        std::vector<boost::shared_ptr<SimpleStimulus> > stimuli_applied_2D;
        HeartConfig::Instance()->GetStimuli(stimuli_applied_2D, stimulated_areas_2D);

        TS_ASSERT_EQUALS(stimuli_applied_2D.size(), 2u);
        TS_ASSERT_EQUALS(stimulated_areas_2D.size(), 2u);
        
        //covering the 1D case
        std::vector<ChasteCuboid<1> > stimulated_areas_1D;
        std::vector<boost::shared_ptr<SimpleStimulus> > stimuli_applied_1D;
        HeartConfig::Instance()->GetStimuli(stimuli_applied_1D, stimulated_areas_1D);

        TS_ASSERT_EQUALS(stimuli_applied_1D.size(), 2u);
        TS_ASSERT_EQUALS(stimulated_areas_1D.size(), 2u);

        TS_ASSERT(stimulated_areas[1].DoesContain(ChastePoint<3>(-2, 0, -2)));
        TS_ASSERT( ! stimulated_areas[1].DoesContain(ChastePoint<3>(-6, -6, -6)));


        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "ChasteResults");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputFilenamePrefix(), "SimulationResults");

        c_vector<double, 3> intra_conductivities;
        HeartConfig::Instance()->GetIntracellularConductivities(intra_conductivities);
        TS_ASSERT_EQUALS(intra_conductivities[0], 1.75);
        TS_ASSERT_EQUALS(intra_conductivities[1], 1.75);
        TS_ASSERT_EQUALS(intra_conductivities[2], 1.75);

         TS_ASSERT_EQUALS(HeartConfig::Instance()->GetBathConductivity(), 7.0);

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

        TS_ASSERT(HeartConfig::Instance()->IsPostProcessingSectionPresent());

        TS_ASSERT(HeartConfig::Instance()->IsApdMapsRequested());
        std::vector<std::pair<double,double> > apd_maps_requested;
        HeartConfig::Instance()->GetApdMaps(apd_maps_requested);
        TS_ASSERT_EQUALS(apd_maps_requested.size(), 1u);
        TS_ASSERT_EQUALS(apd_maps_requested[0].first, 90.0);
        TS_ASSERT_EQUALS(apd_maps_requested[0].second, -30.0);

        TS_ASSERT(HeartConfig::Instance()->IsUpstrokeTimeMapsRequested());
        std::vector<double> upstroke_time_maps_requested;
        HeartConfig::Instance()->GetUpstrokeTimeMaps(upstroke_time_maps_requested);
        TS_ASSERT_EQUALS(upstroke_time_maps_requested.size(), 1u);
        TS_ASSERT_EQUALS(upstroke_time_maps_requested[0], -30.0);

        TS_ASSERT(HeartConfig::Instance()->IsMaxUpstrokeVelocityMapRequested());
        std::vector<double> upstroke_velocity_maps_requested;
        HeartConfig::Instance()->GetMaxUpstrokeVelocityMaps(upstroke_velocity_maps_requested);
        TS_ASSERT_EQUALS(upstroke_velocity_maps_requested.size(), 1u);
        TS_ASSERT_EQUALS(upstroke_velocity_maps_requested[0], -30.0);

        TS_ASSERT(HeartConfig::Instance()->IsConductionVelocityMapsRequested());
        std::vector<unsigned> conduction_velocity_maps_requested;
        HeartConfig::Instance()->GetConductionVelocityMaps(conduction_velocity_maps_requested);
        TS_ASSERT_EQUALS(conduction_velocity_maps_requested.size(), 2u);
        TS_ASSERT_EQUALS(conduction_velocity_maps_requested[0], 10u);
        TS_ASSERT_EQUALS(conduction_velocity_maps_requested[1], 20u);


        /// \todo: refactor from here until the end of the test into a different test
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersLoadMesh.xml");

        TS_ASSERT(!HeartConfig::Instance()->GetCreateMesh());
        TS_ASSERT(HeartConfig::Instance()->GetLoadMesh());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshName(), "foo");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetConductivityMedia(), cp::media_type::NoFibreOrientation);

        //Try reading through an empty parameters file into a fully-populated default
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteEmpty.xml");
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteParametersFullFormat.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
    }
    
    void TestGetHeterogeneities()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");
        ///////////////
        //ionic models
        //////////////
        std::vector<ChasteCuboid<3> > ionic_model_regions;
        std::vector<cp::ionic_model_selection_type> ionic_models_defined;
        HeartConfig::Instance()->GetIonicModelRegions(ionic_model_regions,
                                                      ionic_models_defined);

        TS_ASSERT_EQUALS(ionic_model_regions.size(), 2u);
        TS_ASSERT_EQUALS(ionic_models_defined.size(), 2u);

        TS_ASSERT(ionic_model_regions[0].DoesContain(ChastePoint<3>(-1.95, 0, 0)));
        std::string model_zero("heart/dynamic/libDynamicallyLoadableLr91.so");
        TS_ASSERT(ionic_models_defined[0].Dynamic().present());
        TS_ASSERT_EQUALS(ionic_models_defined[0].Dynamic().get().Path().relative_to(), cp::relative_to_type::chaste_source_root);
        TS_ASSERT_EQUALS(ionic_models_defined[0].Dynamic().get().Path(), model_zero);
        TS_ASSERT(ionic_models_defined[1].Hardcoded().present());
        TS_ASSERT_EQUALS(ionic_models_defined[1].Hardcoded().get(), cp::ionic_models_available_type::DifrancescoNoble);

        //cover the 2D case
        std::vector<ChasteCuboid<2> > ionic_model_regions_2D;
        std::vector<cp::ionic_model_selection_type> ionic_models_defined_2D;
        HeartConfig::Instance()->GetIonicModelRegions(ionic_model_regions_2D,
                                                      ionic_models_defined_2D);

        TS_ASSERT_EQUALS(ionic_model_regions_2D.size(), 2u);
        TS_ASSERT_EQUALS(ionic_models_defined_2D.size(), 2u);

        TS_ASSERT(ionic_model_regions_2D[0].DoesContain(ChastePoint<2>(-1.95, 0)));
        TS_ASSERT(ionic_models_defined_2D[0].Dynamic().present());
        TS_ASSERT_EQUALS(ionic_models_defined_2D[0].Dynamic().get().Path().relative_to(), cp::relative_to_type::chaste_source_root);
        TS_ASSERT_EQUALS(ionic_models_defined_2D[0].Dynamic().get().Path(), model_zero);
        TS_ASSERT(ionic_models_defined_2D[1].Hardcoded().present());
        TS_ASSERT_EQUALS(ionic_models_defined_2D[1].Hardcoded().get(), cp::ionic_models_available_type::DifrancescoNoble);
        
        //cover the 1D case
        std::vector<ChasteCuboid<1> > ionic_model_regions_1D;
        std::vector<cp::ionic_model_selection_type> ionic_models_defined_1D;
        HeartConfig::Instance()->GetIonicModelRegions(ionic_model_regions_1D,
                                                      ionic_models_defined_1D);

        TS_ASSERT_EQUALS(ionic_model_regions_1D.size(), 2u);
        TS_ASSERT_EQUALS(ionic_models_defined_1D.size(), 2u);

        TS_ASSERT(ionic_model_regions_1D[0].DoesContain(ChastePoint<1>(-1.95)));
        TS_ASSERT(ionic_models_defined_1D[0].Dynamic().present());
        TS_ASSERT_EQUALS(ionic_models_defined_1D[0].Dynamic().get().Path().relative_to(), cp::relative_to_type::chaste_source_root);
        TS_ASSERT_EQUALS(ionic_models_defined_1D[0].Dynamic().get().Path(), model_zero);
        TS_ASSERT(ionic_models_defined_1D[1].Hardcoded().present());
        TS_ASSERT_EQUALS(ionic_models_defined_1D[1].Hardcoded().get(), cp::ionic_models_available_type::DifrancescoNoble);
        
        ///////////////
        //Cell heterogeneities
        //////////////
        
        std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas_3D;
        std::vector<double> scale_factor_gks_3D;
        std::vector<double> scale_factor_ito_3D;
        std::vector<double> scale_factor_gkr_3D;  
        
        TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularTransmuralHeterogeneitiesRequested(), false);
        
        HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas_3D,
                                                        scale_factor_gks_3D,
                                                        scale_factor_ito_3D,
                                                        scale_factor_gkr_3D);
                                                        
        TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularlHeterogeneitiesSpecifiedByCuboids(), true);
        TS_ASSERT(cell_heterogeneity_areas_3D[0]->DoesContain(ChastePoint<3>(-1.0, 0, 0)));
        TS_ASSERT_EQUALS(scale_factor_gks_3D[0], 0.462);
        TS_ASSERT_EQUALS(scale_factor_ito_3D[0], 0.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_3D[0], 1.0);
        TS_ASSERT_EQUALS(scale_factor_gks_3D[1], 1.154);
        TS_ASSERT_EQUALS(scale_factor_ito_3D[1], 0.85);
        TS_ASSERT_EQUALS(scale_factor_gkr_3D[1], 1.0);
        
        for (unsigned i = 0; i<cell_heterogeneity_areas_3D.size();i++)
        {
            delete cell_heterogeneity_areas_3D[i];
        }

        //cover the 2D case 
        std::vector<AbstractChasteRegion<2>* > cell_heterogeneity_areas_2D;
        std::vector<double> scale_factor_gks_2D;
        std::vector<double> scale_factor_ito_2D;
        std::vector<double> scale_factor_gkr_2D;
        HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas_2D,
                                                        scale_factor_gks_2D,
                                                        scale_factor_ito_2D,
                                                        scale_factor_gkr_2D);

        TS_ASSERT(cell_heterogeneity_areas_2D[0]->DoesContain(ChastePoint<2>(-1.0, 0)));
        TS_ASSERT_EQUALS(scale_factor_gks_2D[0], 0.462);
        TS_ASSERT_EQUALS(scale_factor_ito_2D[0], 0.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_2D[0], 1.0);
        TS_ASSERT_EQUALS(scale_factor_gks_2D[1], 1.154);
        TS_ASSERT_EQUALS(scale_factor_ito_2D[1], 0.85);
        TS_ASSERT_EQUALS(scale_factor_gkr_2D[1], 1.0);
        
        for (unsigned i = 0; i<cell_heterogeneity_areas_2D.size();i++)
        {
            delete cell_heterogeneity_areas_2D[i];
        }
        
         //cover the 1D case 
        std::vector<AbstractChasteRegion<1>* > cell_heterogeneity_areas_1D;
        std::vector<double> scale_factor_gks_1D;
        std::vector<double> scale_factor_ito_1D;
        std::vector<double> scale_factor_gkr_1D;
        HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas_1D,
                                                        scale_factor_gks_1D,
                                                        scale_factor_ito_1D,
                                                        scale_factor_gkr_1D);

        TS_ASSERT(cell_heterogeneity_areas_1D[0]->DoesContain(ChastePoint<1>(-1.0)));
        TS_ASSERT_EQUALS(scale_factor_gks_1D[0], 0.462);
        TS_ASSERT_EQUALS(scale_factor_ito_1D[0], 0.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_1D[0], 1.0);
        TS_ASSERT_EQUALS(scale_factor_gks_1D[1], 1.154);
        TS_ASSERT_EQUALS(scale_factor_ito_1D[1], 0.85);
        TS_ASSERT_EQUALS(scale_factor_gkr_1D[1], 1.0);
        
        for (unsigned i = 0; i<cell_heterogeneity_areas_1D.size();i++)
        {
            delete cell_heterogeneity_areas_1D[i];
        }
   
        ///////////////
        //Conductivity heterogeneities
        //////////////     
        
        std::vector<ChasteCuboid<3> > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        TS_ASSERT(conductivities_heterogeneity_areas[0].DoesContain(ChastePoint<3>(1.95, 0, 0)));
        TS_ASSERT_EQUALS(intra_h_conductivities[0][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities[0][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities[1][0], 0.75);
        
        //cover the 2D case
        std::vector<ChasteCuboid<2> > conductivities_heterogeneity_areas_2D;
        std::vector< c_vector<double,3> > intra_h_conductivities_2D;
        std::vector< c_vector<double,3> > extra_h_conductivities_2D;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_2D,
                                                                intra_h_conductivities_2D,
                                                                extra_h_conductivities_2D);

        TS_ASSERT(conductivities_heterogeneity_areas_2D[0].DoesContain(ChastePoint<2>(1.95, 0)));
        TS_ASSERT_EQUALS(intra_h_conductivities_2D[0][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities_2D[0][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities_2D[1][0], 0.75);
        
        //cover the 1D case
        std::vector<ChasteCuboid<1> > conductivities_heterogeneity_areas_1D;
        std::vector< c_vector<double,3> > intra_h_conductivities_1D;
        std::vector< c_vector<double,3> > extra_h_conductivities_1D;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_1D,
                                                                intra_h_conductivities_1D,
                                                                extra_h_conductivities_1D);

        TS_ASSERT(conductivities_heterogeneity_areas_1D[0].DoesContain(ChastePoint<1>(1.95)));
        TS_ASSERT_EQUALS(intra_h_conductivities_1D[0][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities_1D[0][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities_1D[1][0], 0.75);    
        
        //extracellular conductivities
        c_vector<double, 3> extra_conductivities;
        HeartConfig::Instance()->GetExtracellularConductivities(extra_conductivities);
        TS_ASSERT_EQUALS(extra_h_conductivities[1][0], extra_conductivities[0]); 
        
        TS_ASSERT_EQUALS(extra_conductivities[0], 7.0);
        TS_ASSERT_EQUALS(extra_conductivities[1], 7.0);
        TS_ASSERT_EQUALS(extra_conductivities[2], 7.0);  
    }

    void TestIsMeshProvided()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteEmpty.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsMeshProvided(), false);
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetLoadMesh(),"No Mesh provided (neither default nor user defined)");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCreateMesh(), "No Mesh provided (neither default nor user defined)");

        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersLoadMesh.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsMeshProvided(), true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetLoadMesh(), true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCreateMesh(), false);
    }
    
    void TestTransmuralHeterogeneities()
    {
        {
            HeartConfig::Reset();
            //the _unsupported file has valid transmural heterogeneity definition for cellular heterogeneities, but transmural heterogeneities defined for other things we don't support yet.
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersCellHeterogeneities_unsupported.xml");
            
            std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
            std::vector<double> scale_factor_gks;
            std::vector<double> scale_factor_ito;
            std::vector<double> scale_factor_gkr;
            
            //before we call GetCellHeterogeneities, the flag should be false
            TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularTransmuralHeterogeneitiesRequested(), false);
            
            //and the indices with their initial values
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEpiLayerIndex(), UINT_MAX-3u);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMidLayerIndex(), UINT_MAX-3u);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEndoLayerIndex(), UINT_MAX-3u);
            
            HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                            scale_factor_gks,
                                                            scale_factor_ito,
                                                            scale_factor_gkr);
                                                            
            //in this file, they are supplied as epi first, then endo, then mid
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEpiLayerIndex(), 0u);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMidLayerIndex(), 2u);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEndoLayerIndex(), 1u);
            
                                                                       
            //now the flag for transmural cellular heterogeneities should be true
            TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularTransmuralHeterogeneitiesRequested(), true);
                                                          
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEpiLayerFraction(),0.3);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEndoLayerFraction(),0.3);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMidLayerFraction(),0.4);
            TS_ASSERT_EQUALS(scale_factor_gks[0], 0.462);
            TS_ASSERT_EQUALS(scale_factor_ito[0], 0.0);
            TS_ASSERT_EQUALS(scale_factor_gkr[0], 1.0);
            
            //covering the case when the user specify transmural layers for conductivities (not supported and probably not worth considering)... 
            std::vector<ChasteCuboid<3> > conductivities_heterogeneity_areas;
            std::vector< c_vector<double,3> > intra_h_conductivities;
            std::vector< c_vector<double,3> > extra_h_conductivities;
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                                          intra_h_conductivities,
                                                                                          extra_h_conductivities),
                                  "Definition of transmural layers is not allowed for conductivities heterogeneities, you may use fibre orientation support instead");
             
             //covering the case when the user specify transmural layers for stimulated areas (not yet supported)... 
            std::vector<boost::shared_ptr<SimpleStimulus> > stimuli_applied;
            std::vector<ChasteCuboid<3> > stimulated_area;
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_area),
                                  "Definition of transmural layers is not yet supported for specifying stimulated areas, please use cuboids instead");
            
            //covering the case when the user specify transmural layers for ionic model heterogeneities (not yet supported)... 
            std::vector<ChasteCuboid<3> > ionic_model_regions;
            std::vector<cp::ionic_model_selection_type> ionic_models;
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetIonicModelRegions(ionic_model_regions, ionic_models),
                                  "Definition of transmural layers is not yet supported for defining different ionic models, please use cuboids instead");
        }
        //covers the case when the user supplies numbers that do not add up to 1
        {
            HeartConfig::Reset();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersCellHeterogeneities_inconsistent.xml");
            
            std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
            std::vector<double> scale_factor_gks;
            std::vector<double> scale_factor_ito;
            std::vector<double> scale_factor_gkr;
            
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                                                  scale_factor_gks,
                                                                                  scale_factor_ito,
                                                                                  scale_factor_gkr),
                                  "Summation of epicardial, midmyocardial and endocardial fractions should be 1");
                                                            
        }
        //covers the case when the user supplies negative numbers
        {
            HeartConfig::Reset();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersCellHeterogeneities_negative.xml");
            
            std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
            std::vector<double> scale_factor_gks;
            std::vector<double> scale_factor_ito;
            std::vector<double> scale_factor_gkr;
            
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                                                  scale_factor_gks,
                                                                                  scale_factor_ito,
                                                                                  scale_factor_gkr),
                                  "Fractions must be positive");
                                                            
        }
        //covers the case when the user supplies only two layers
        {
            HeartConfig::Reset();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersCellHeterogeneities_incomplete.xml");
            
            std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
            std::vector<double> scale_factor_gks;
            std::vector<double> scale_factor_ito;
            std::vector<double> scale_factor_gkr;

            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                                                  scale_factor_gks,
                                                                                  scale_factor_ito,
                                                                                  scale_factor_gkr),
                                  "Three specifications of layers must be supplied");
            
 
            TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularTransmuralHeterogeneitiesRequested(), true);           
            //only epi and endo, in this order, are supplied in this file
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEpiLayerIndex(), 0u);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMidLayerIndex(), UINT_MAX-3u);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetEndoLayerIndex(), 1u);           
                                                            
        }
        //cuboids and layers together are not yet supported
        {
            HeartConfig::Reset();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersCellHeterogeneities_cuboids_and_layers.xml");
            
            std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
            std::vector<double> scale_factor_gks;
            std::vector<double> scale_factor_ito;
            std::vector<double> scale_factor_gkr;
            
            TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularlHeterogeneitiesSpecifiedByCuboids(), false);        
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                                                  scale_factor_gks,
                                                                                  scale_factor_ito,
                                                                                  scale_factor_gkr),
                                  "Specification of cellular heterogeneities by cuboids and layers at the same time is not yet supported");
                                  
            TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularlHeterogeneitiesSpecifiedByCuboids(), true);           
        }
        
    }
    void Test2dProblems()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParameters2D.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSpaceDimension(),
                         2u);

        TS_ASSERT(HeartConfig::Instance()->GetCreateMesh());
        TS_ASSERT(!HeartConfig::Instance()->GetLoadMesh());

        c_vector<double, 3> slab_dimensions;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetSlabDimensions(slab_dimensions),
                "Tissue slabs can only be defined in 3D");
        c_vector<double, 1> fibre_length;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetFibreLength(fibre_length),
                "Tissue fibres can only be defined in 1D");

        c_vector<double, 2> sheet_dimensions;
        HeartConfig::Instance()->GetSheetDimensions(sheet_dimensions);

        TS_ASSERT_EQUALS(sheet_dimensions[0], 4.0);
        TS_ASSERT_EQUALS(sheet_dimensions[1], 0.1);

        double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();
        TS_ASSERT_EQUALS(inter_node_space, 0.1);
    }

    void Test1dProblems()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParameters1D.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSpaceDimension(),
                         1u);

        TS_ASSERT(HeartConfig::Instance()->GetCreateMesh());
        TS_ASSERT(!HeartConfig::Instance()->GetLoadMesh());

        c_vector<double, 3> slab_dimensions;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetSlabDimensions(slab_dimensions),
                "Tissue slabs can only be defined in 3D");
        c_vector<double, 2> sheet_dimensions;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetSheetDimensions(sheet_dimensions),
                "Tissue sheets can only be defined in 2D");

        c_vector<double, 1> fibre_length;
        HeartConfig::Instance()->GetFibreLength(fibre_length);

        TS_ASSERT_EQUALS(fibre_length[0], 4.0);

        double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();
        TS_ASSERT_EQUALS(inter_node_space, 0.1);
    }


    void TestSetFunctions() throw(Exception)
    {
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingSectionPresent(), true);
        HeartConfig::Instance()->SetSimulationDuration(35.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 35.0);

        HeartConfig::Instance()->SetDomain(cp::domain_type::Bi);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDomain(), cp::domain_type::Bi);

        HeartConfig::Instance()->SetDefaultIonicModel(cp::ionic_models_available_type::LuoRudyIBackwardEuler);
        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(), cp::ionic_models_available_type::LuoRudyIBackwardEuler);

        HeartConfig::Instance()->SetDefaultIonicModel(cp::ionic_models_available_type::MahajanShiferaw);
        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(), cp::ionic_models_available_type::MahajanShiferaw);

        HeartConfig::Instance()->SetDefaultIonicModel(cp::ionic_models_available_type::HodgkinHuxley);
        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(), cp::ionic_models_available_type::HodgkinHuxley);

        HeartConfig::Instance()->SetDefaultIonicModel(cp::ionic_models_available_type::tenTusscher2006);
        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(), cp::ionic_models_available_type::tenTusscher2006);

        HeartConfig::Instance()->SetDefaultIonicModel(cp::ionic_models_available_type::DifrancescoNoble);
        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(), cp::ionic_models_available_type::DifrancescoNoble);

        TS_ASSERT(!HeartConfig::Instance()->GetConductivityHeterogeneitiesProvided());

        std::vector< c_vector<double,3> > intraConductivities;
        std::vector< c_vector<double,3> > extraConductivities;

        std::vector<ChasteCuboid<3> > input_areas;
        ChastePoint<3> lower1(-1.0, -1.0, -1.0);
        ChastePoint<3> upper1( 1.0,  1.0,  1.0);
        input_areas.push_back(ChasteCuboid<3> (lower1, upper1));
        intraConductivities.push_back( Create_c_vector(2.5, 2.5, 2.5) );
        extraConductivities.push_back( Create_c_vector(8.5, 8.5, 8.5) );

        ChastePoint<3> lower2(-2.0, -2.0, -2.0);
        ChastePoint<3> upper2(-1.0, -1.0, -1.0);
        input_areas.push_back(ChasteCuboid<3> (lower2, upper2));
        intraConductivities.push_back( Create_c_vector(1.0, 0.5, 0.4) );
        extraConductivities.push_back( Create_c_vector(7.0, 6.5, 6.4) );

        HeartConfig::Instance()->SetConductivityHeterogeneities(input_areas, intraConductivities, extraConductivities);

        TS_ASSERT(HeartConfig::Instance()->GetConductivityHeterogeneitiesProvided());

        std::vector<ChasteCuboid<3> > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas.size(), 2u);
        TS_ASSERT_EQUALS(intra_h_conductivities.size(), 2u);
        TS_ASSERT_EQUALS(extra_h_conductivities.size(), 2u);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas[0].DoesContain(ChastePoint<3>(0.0, 0.0, 0.0)), true);
        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas[0].DoesContain(ChastePoint<3>(-1.5, -1.5, -1.5)), false);
        TS_ASSERT_EQUALS(intra_h_conductivities[0][0], 2.5);
        TS_ASSERT_EQUALS(extra_h_conductivities[0][0], 8.5);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas[1].DoesContain(ChastePoint<3>(-1.5, -1.5, -1.5)), true);
        TS_ASSERT_EQUALS(intra_h_conductivities[1][0], 1.0);


        std::vector<ChasteCuboid<3> > ionic_model_regions;
        std::vector<cp::ionic_model_selection_type> ionic_models;

        //No ionic model regions
        HeartConfig::Instance()->GetIonicModelRegions(ionic_model_regions,
                                                      ionic_models);
        TS_ASSERT_EQUALS(ionic_model_regions.size(), 0u);

        //Set ionic regions
        std::vector<cp::ionic_model_selection_type> input_ionic_models;
        cp::ionic_model_selection_type model1;
        cp::ionic_models_available_type model1_type = cp::ionic_models_available_type::HodgkinHuxley;
        model1.Hardcoded(model1_type);
        input_ionic_models.push_back(model1);
        cp::ionic_model_selection_type model2;
        std::string model2_path_str("made-up-path");
        cp::path_type model2_path(model2_path_str);
        model2.Dynamic(model2_path);
        input_ionic_models.push_back(model2);
        HeartConfig::Instance()->SetIonicModelRegions(input_areas, input_ionic_models);

        //Now there are 2 ionic model regions
        HeartConfig::Instance()->GetIonicModelRegions(ionic_model_regions,
                                                      ionic_models);
        TS_ASSERT_EQUALS(ionic_model_regions.size(), 2u);
        // Check they were set correctly
        TS_ASSERT(ionic_models[0].Hardcoded().present());
        TS_ASSERT_EQUALS(ionic_models[0].Hardcoded().get(), model1_type);
        TS_ASSERT(ionic_models[1].Dynamic().present());
        TS_ASSERT_EQUALS(ionic_models[1].Dynamic().get().Path(), model2_path_str);
        TS_ASSERT_EQUALS(ionic_models[1].Dynamic().get().Path().relative_to(), cp::relative_to_type::cwd);

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

        //One-dimensional set
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(-3.0));
        c_vector<double, 3> extra;
        HeartConfig::Instance()->GetExtracellularConductivities(extra);
        TS_ASSERT_EQUALS(extra[0], -3.0);
        TS_ASSERT_EQUALS(extra[1], 0.0);
        TS_ASSERT_EQUALS(extra[2], 0.0);

        //Two-dimensional set
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(-3.0, -2.0));
        HeartConfig::Instance()->GetExtracellularConductivities(extra);
        TS_ASSERT_EQUALS(extra[0], -3.0);
        TS_ASSERT_EQUALS(extra[1], -2.0);
        TS_ASSERT_EQUALS(extra[2], 0.0);

        //Three-dimensional set
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(-3.0, -2.0, -1.0));
        HeartConfig::Instance()->GetExtracellularConductivities(extra);
        TS_ASSERT_EQUALS(extra[0], -3.0);
        TS_ASSERT_EQUALS(extra[1], -2.0);
        TS_ASSERT_EQUALS(extra[2], -1.0);

        HeartConfig::Instance()->SetBathConductivity(150);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetBathConductivity(), 150);

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
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdeTimeStep(0.2), "Ode time-step should not be greater than pde time-step");

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(-0.1, 0.1, 0.1), "Ode time-step should be positive");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, -0.1, 0.1), "Pde time-step should be positive");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, -0.1), "Printing time-step should be positive");

        //Throws when we try to print more often than the pde time step
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.1), "Printing time-step should not be smaller than PDE time step");
         //Throws when printing step is not a multiple of pde time step
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.3), "Printing time-step should be a multiple of PDE time step");

        // Throws because ode time step is bigger than pde time step
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.3), "Printing time-step should be a multiple of PDE time step");

        /*
         *  Set ODE, PDE, and printing timestep to something meaningful and test SetCheckpointSimulation() exceptions. 
         */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.4);
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetCheckpointSimulation(true, 1.0, 3), "Checkpoint time-step should be a multiple of printing time step");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetCheckpointSimulation(true, -2.0, 3), "Printing time-step should be positive");

        HeartConfig::Instance()->SetUseRelativeTolerance(1e-4);
        TS_ASSERT(HeartConfig::Instance()->GetUseRelativeTolerance());
        TS_ASSERT(HeartConfig::Instance()->GetUseAbsoluteTolerance() == false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetRelativeTolerance(), 1e-4);

        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-11);
        TS_ASSERT(HeartConfig::Instance()->GetUseAbsoluteTolerance());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetAbsoluteTolerance(), 1e-11);
        //Check that relative tolerance is disabled
        TS_ASSERT(HeartConfig::Instance()->GetUseRelativeTolerance() == false);
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetRelativeTolerance(),
                "Relative tolerance is not set in Chaste parameters");


        HeartConfig::Instance()->SetUseRelativeTolerance(1e-4);
        TS_ASSERT(HeartConfig::Instance()->GetUseRelativeTolerance());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetRelativeTolerance(), 1e-4);
        //Check that absolute tolerance is disabled
        TS_ASSERT(HeartConfig::Instance()->GetUseAbsoluteTolerance() == false);
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetAbsoluteTolerance(),
                "Absolute tolerance is not set in Chaste parameters");

        HeartConfig::Instance()->SetKSPSolver("cg");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "cg")==0);

        HeartConfig::Instance()->SetKSPSolver("gmres");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "gmres")==0);

        HeartConfig::Instance()->SetKSPSolver("symmlq");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "symmlq")==0);

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetKSPSolver("foobar"),"Unknown solver type provided");

        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "jacobi")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "bjacobi")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("ilu");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "ilu")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("hypre");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "hypre")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("ml");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "ml")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("spai");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "spai")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("blockdiagonal");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "blockdiagonal")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("ldufactorisation");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "ldufactorisation")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("none");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "none")==0);

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetKSPPreconditioner("foobar"),
                "Unknown preconditioner type provided");

        // Tests for set functions of postprocessing

        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsApdMapsRequested(), false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingRequested(), false);
        std::vector<std::pair<double,double> > apds, apd_maps;
        apds.push_back(std::pair<double, double>(90,-30));//reploarisation percentage first, as per schema
        HeartConfig::Instance()->SetApdMaps(apds);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsApdMapsRequested(), true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingRequested(), true);
        HeartConfig::Instance()->GetApdMaps(apd_maps);
        TS_ASSERT_EQUALS(apd_maps.size(),1u);
        TS_ASSERT_EQUALS(apd_maps[0].first,90);
        TS_ASSERT_EQUALS(apd_maps[0].second,-30);

        apds[0].first = 80;//reploarisation percentage first, as per schema
        apds[0].second = -45;
        HeartConfig::Instance()->SetApdMaps(apds);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsApdMapsRequested(), true);

        HeartConfig::Instance()->GetApdMaps(apd_maps);
        TS_ASSERT_EQUALS(apd_maps.size(),1u);
        TS_ASSERT_EQUALS(apd_maps[0].first,80);
        TS_ASSERT_EQUALS(apd_maps[0].second,-45);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsUpstrokeTimeMapsRequested(), false);
        std::vector<double> upstroke_time_map, upstroke_time_map_get;
        upstroke_time_map.push_back(25.0);
        upstroke_time_map.push_back(55.0);
        HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsUpstrokeTimeMapsRequested(), true);
        HeartConfig::Instance()->GetUpstrokeTimeMaps(upstroke_time_map_get);
        TS_ASSERT_EQUALS(upstroke_time_map_get.size(),2u);
        TS_ASSERT_EQUALS(upstroke_time_map_get[0],25);
        TS_ASSERT_EQUALS(upstroke_time_map_get[1],55);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsMaxUpstrokeVelocityMapRequested(), false);
        std::vector<double> upstroke_velocity_map, upstroke_velocity_map_get;
        upstroke_velocity_map.push_back(25.0);
        upstroke_velocity_map.push_back(55.0);
        HeartConfig::Instance()->SetMaxUpstrokeVelocityMaps(upstroke_velocity_map);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsMaxUpstrokeVelocityMapRequested(), true);
        HeartConfig::Instance()->GetMaxUpstrokeVelocityMaps(upstroke_velocity_map_get);
        TS_ASSERT_EQUALS(upstroke_velocity_map_get.size(),2u);
        TS_ASSERT_EQUALS(upstroke_velocity_map_get[0],25);
        TS_ASSERT_EQUALS(upstroke_velocity_map_get[1],55);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsConductionVelocityMapsRequested(), false);
        std::vector<unsigned> conduction_velocity_map, conduction_velocity_map_get;
        conduction_velocity_map.push_back(25u);
        HeartConfig::Instance()->SetConductionVelocityMaps(conduction_velocity_map);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsConductionVelocityMapsRequested(), true);
        HeartConfig::Instance()->GetConductionVelocityMaps(conduction_velocity_map_get);
        TS_ASSERT_EQUALS(conduction_velocity_map_get.size(),1u);
        TS_ASSERT_EQUALS(conduction_velocity_map_get[0],25u);
    }

    void TestWrite() throw(Exception)
    {
        OutputFileHandler output_file_handler("Xml/output", true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 0.01);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(1.1,2.2,4.4);
        HeartConfig::Instance()->SetOutputDirectory("Xml");
        HeartConfig::Instance()->Write();

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 1.1);
        HeartConfig::Reset();
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 0.01);
        //Reload the other XML
        HeartConfig::Instance()->SetParametersFile(output_file_handler.GetOutputDirectoryFullPath()+"ChasteParameters.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOdeTimeStep(), 1.1);
    }

    void TestArchiving()
    {
        //Archive
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        handler.SetArchiveDirectory();
        archive_filename = handler.GetOutputDirectoryFullPath() + "heart_config.arch";

        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");

        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        cp::ionic_models_available_type user_ionic = HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get();
        TS_ASSERT( user_ionic == cp::ionic_models_available_type::FaberRudy2000 );
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);

        std::ofstream ofs(archive_filename.c_str());
        boost::archive::text_oarchive output_arch(ofs);
        HeartConfig* const p_archive_heart_config = HeartConfig::Instance();
        output_arch << static_cast<const HeartConfig&>(*p_archive_heart_config);

        ofs.close();

        HeartConfig::Reset();

        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(), cp::ionic_models_available_type::LuoRudyI);

        // We split the two load attempts into their own scopes to avoid a
        // memory leak (uninitialised value).
        {
            // Try a load that fails
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            HeartConfig* p_heart_config = HeartConfig::Instance();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersResumeSimulationWrongDimension.xml");
            TS_ASSERT_THROWS_THIS(input_arch >> (*p_heart_config), "Problem type and space dimension should match when restarting a simulation.");
        }

        {
            // Try a correct load
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            HeartConfig* p_heart_config = HeartConfig::Instance();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersResumeSimulation.xml");
            input_arch >> (*p_heart_config);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 20.0);
            TS_ASSERT(p_heart_config->GetDefaultIonicModel().Hardcoded().present());
            TS_ASSERT_EQUALS( user_ionic, p_heart_config->GetDefaultIonicModel().Hardcoded().get());
        }
        
        {
            // Try a correct load without checkpointing to make sure it's not inherited 
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            HeartConfig* p_heart_config = HeartConfig::Instance();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersResumeSimulationWithoutCheckpointing.xml");
            input_arch >> (*p_heart_config);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 20.0);
            TS_ASSERT(p_heart_config->GetDefaultIonicModel().Hardcoded().present());
            TS_ASSERT_EQUALS( user_ionic, p_heart_config->GetDefaultIonicModel().Hardcoded().get());
            TS_ASSERT(!HeartConfig::Instance()->GetCheckpointSimulation());
        }
    }

    void TestExceptions() throw (Exception)
    {
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetDefaultsFile("DoesNotExist.xml"),
                "Missing file parsing configuration file: DoesNotExist.xml");
        HeartConfig::Reset();
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteInconsistent.xml"),
                "Ode time-step should not be greater than pde time-step");
        HeartConfig::Reset();
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteInconsistent.xml"),
                "Ode time-step should not be greater than pde time-step");
        HeartConfig::Reset();
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteInconsistentCheckpointTimestep.xml"),
                "Checkpoint time-step should be a multiple of printing time step");

        //Can't open a directory
        HeartConfig::Instance()->SetOutputDirectory("../../../");
        TS_ASSERT_THROWS_CONTAINS(HeartConfig::Instance()->Write(), "due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");

        // Can't open a file for writing
        std::string command = OutputFileHandler::GetChasteTestOutputDirectory() + "no_write_access";
        mkdir(command.c_str(), 0444);
        HeartConfig::Instance()->SetOutputDirectory("no_write_access");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->Write(),"Could not open XML file in HeartConfig");
        chmod(command.c_str(), 0755);
        rmdir(command.c_str());

        // A non-parsing exception, for coverage. Very hard to get in practice!
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/EmptyRoot.xml"),
                              "XML parsing error in configuration file: heart/test/data/xml/EmptyRoot.xml");
    }

    /**
     *  The following test is aimed at checking that the ChasteParameters.xml file,
     *  which is distributed with the executable, remains valid.
     */
    void TestChasteParametersFile() throw (Exception)
    {
        HeartConfig::Instance()->SetParametersFile("ChasteParameters.xml");
    }

    /**
     * And here we try to check that using old XML or XSD files does The Right Thing.
     */
    void TestVersioning()
    {
        // Broken schema should throw
        HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/BrokenSchema.xml"),
                "XML parsing error in configuration file: heart/test/data/xml/BrokenSchema.xml");
        // But if we use the fixed schema location, it should be OK.
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        TS_ASSERT_THROWS_NOTHING(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/BrokenSchema.xml"));

        // Can release 1 xml be loaded with release 1 schema?
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease1.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease1.xml");

        // Check that release 1 xml can be loaded with release 1.1 schema
        HeartConfig::Reset();
        HeartConfig::SchemaLocationsMap schema_locations;
        schema_locations[""] = std::string(ChasteBuildInfo::GetRootDir()) + "/heart/test/data/xml/ChasteParametersRelease1_1.xsd";
        HeartConfig::Instance()->SetFixedSchemaLocations(schema_locations);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease1.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease1.xml");

        // Check that release 1 xml can be loaded with latest schema
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease1.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease1.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingSectionPresent(), false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingRequested(), false);

        // Check that release 1.1 xml can be loaded with release 1.1 schema
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease1_1.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease1_1.xml");

        // Check that release 1.1 xml can be loaded with latest schema
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease1_1.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease1_1.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingSectionPresent(), true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingRequested(), false);
    }

    /**
     * Test whether we can use a schema that has spaces in the pathname.
     * This gives some indication of whether Chaste will cope being checked out into
     * a path with spaces.
     */
    void TestSpacesInPath()
    {
        HeartConfig::Reset();
        HeartConfig::SchemaLocationsMap schema_locations;
        schema_locations[""] = std::string(ChasteBuildInfo::GetRootDir()) + "/heart/test/data/xml/schema with spaces.xsd";
        HeartConfig::Instance()->SetFixedSchemaLocations(schema_locations);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease1_1.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease1_1.xml");
    }

    void TestGetOuputVariablesFromXML()
    {
        // Use the configuration file we just modified.
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");

        // We want a method to check if the user is interested in any extra variable
        TS_ASSERT(HeartConfig::Instance()->GetOutputVariablesProvided());

        // Get them
        std::vector<std::string> output_variables;
        HeartConfig::Instance()->GetOutputVariables(output_variables);

        bool three_variables_defined = (output_variables.size() == 3u);

        // Test three variables were provided
        TS_ASSERT(three_variables_defined);

        // Test the actual names
        if (three_variables_defined)
        {
            TS_ASSERT_EQUALS(output_variables[0],"CaI");
            TS_ASSERT_EQUALS(output_variables[1],"Nai");
            TS_ASSERT_EQUALS(output_variables[2],"Ki");
        }
    }

    void TestSetAndGetOuputVariables()
    {
        // Set the variables we are interested in writing.
        std::vector<std::string> output_variables;
        output_variables.push_back("CaI");
        output_variables.push_back("Nai");
        output_variables.push_back("Ki");

        HeartConfig::Instance()->SetOutputVariables( output_variables );

        // We want a method to check if the user is interested in any extra variable
        TS_ASSERT(HeartConfig::Instance()->GetOutputVariablesProvided());

        // Get them
        std::vector<std::string> got_output_variables;
        HeartConfig::Instance()->GetOutputVariables(got_output_variables);

        bool three_variables_defined = (got_output_variables.size() == 3u);

        // Test three variables were provided
        TS_ASSERT(three_variables_defined);

        // Test the actual names
        if (three_variables_defined)
        {
            TS_ASSERT_EQUALS(got_output_variables[0],"CaI");
            TS_ASSERT_EQUALS(got_output_variables[1],"Nai");
            TS_ASSERT_EQUALS(got_output_variables[2],"Ki");
        }
    }

    void TestSetAndGetArchivingStuff()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");
        TS_ASSERT(HeartConfig::Instance()->IsSimulationDefined());
        TS_ASSERT(!HeartConfig::Instance()->IsSimulationResumed());

        TS_ASSERT(HeartConfig::Instance()->GetCheckpointSimulation());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCheckpointTimestep(),20.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMaxCheckpointsOnDisk(),3u);
        HeartConfig::Instance()->SetCheckpointSimulation(false);
        TS_ASSERT(!HeartConfig::Instance()->GetCheckpointSimulation());
        HeartConfig::Instance()->SetCheckpointSimulation(true, 10.0, 4u);
        TS_ASSERT(HeartConfig::Instance()->GetCheckpointSimulation());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCheckpointTimestep(),10.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMaxCheckpointsOnDisk(),4u);

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetArchivedSimulationDir(),
                              "GetArchivedSimulationDir information is not available in a standard (non-resumed) simulation.");

        // Get the singleton in a clean state
        HeartConfig::Reset();

        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersResumeSimulation.xml");
        TS_ASSERT(!HeartConfig::Instance()->IsSimulationDefined());
        TS_ASSERT(HeartConfig::Instance()->IsSimulationResumed());

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetArchivedSimulationDir(), "ChasteResults_10ms");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 20.0);
        TS_ASSERT(HeartConfig::Instance()->GetCheckpointSimulation());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMaxCheckpointsOnDisk(),3u);        

        // Cover loads of methods where we ask for information that is not present in a ResumedSimulation
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetDefaultIonicModel(), "DefaultIonicModel information is not available in a resumed simulation.")

        std::vector<ChasteCuboid<3> > definedRegions;
        std::vector<cp::ionic_model_selection_type> ionic_models;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetIonicModelRegions(definedRegions,ionic_models),
                              "IonicModelRegions information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->IsMeshProvided(), "Mesh information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCreateMesh(), "Mesh information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCreateSlab(), "Mesh information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCreateSheet(), "Mesh information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCreateFibre(), "Mesh information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetLoadMesh(), "Mesh information is not available in a resumed simulation.");

        c_vector<double, 3> slabDimensions;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetSlabDimensions(slabDimensions),
                              "Slab information is not available in a resumed simulation.");
        c_vector<double, 2> sheet_dimensions;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetSheetDimensions(sheet_dimensions),
                              "Sheet information is not available in a resumed simulation.");
        c_vector<double, 1> fibre_dimensions;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetFibreLength(fibre_dimensions),
                              "Fibre information is not available in a resumed simulation.");

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetInterNodeSpace(), "InterNodeSpace information is not available in a resumed simulation.")
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetMeshName(), "LoadMesh information is not available in a resumed simulation.")
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetConductivityMedia(), "LoadMesh information is not available in a resumed simulation.")

        std::vector<boost::shared_ptr<SimpleStimulus> > stimuli_applied;
        std::vector<ChasteCuboid<3> > stimulated_area;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_area), "Stimuli information is not available in a resumed simulation.")

        std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
        std::vector<double> scale_factor_gks;
        std::vector<double> scale_factor_ito;
        std::vector<double> scale_factor_gkr;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas, scale_factor_gks, scale_factor_ito, scale_factor_gkr),
                              "CellHeterogeneities information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetConductivityHeterogeneitiesProvided(),
                              "ConductivityHeterogeneities information is not available in a resumed simulation.");

        for (unsigned i = 0; i<cell_heterogeneity_areas.size();i++)
        {
            delete cell_heterogeneity_areas[i];
        }

        std::vector<ChasteCuboid<3> > conductivitiesHeterogeneityAreas;
        std::vector< c_vector<double,3> > intraConductivities;
        std::vector< c_vector<double,3> > extraConductivities;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetConductivityHeterogeneities(conductivitiesHeterogeneityAreas,
                                                                                      intraConductivities, extraConductivities),
                              "ConductivityHeterogeneities information is not available in a resumed simulation.");

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetOutputDirectory(),
                              "Simulation/OutputDirectory information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetOutputFilenamePrefix(),
                              "Simulation/OutputFilenamePrefix information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetOutputVariablesProvided(),
                              "OutputVariables information is not available in a resumed simulation.");

        std::vector<std::string> output_variables;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetOutputVariables(output_variables),
                              "OutputVariables information is not available in a resumed simulation.");
    }

    void TestOutputVisualizerSettings()
    {
        // Defaults file doesn't have the OutputVisualizer element
        TS_ASSERT( ! HeartConfig::Instance()->mpDefaultParameters->Simulation().get().OutputVisualizer().present());

        // Parameters file which doesn't specify OutputVisualizer
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteEmpty.xml");

        // Read the user parameters directly - again element missing
        TS_ASSERT( ! HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().present());

        // And the normal Get methods
        TS_ASSERT( HeartConfig::Instance()->GetVisualizeWithMeshalyzer() );
        TS_ASSERT( ! HeartConfig::Instance()->GetVisualizeWithCmgui() );
        TS_ASSERT( ! HeartConfig::Instance()->GetVisualizeWithVtk() );

        // Set methods
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
        TS_ASSERT( ! HeartConfig::Instance()->GetVisualizeWithMeshalyzer() );
        TS_ASSERT_EQUALS(HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().get().meshalyzer(),
                         cp::yesno_type::no);

        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        TS_ASSERT( HeartConfig::Instance()->GetVisualizeWithCmgui() );
        TS_ASSERT_EQUALS(HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().get().cmgui(),
                         cp::yesno_type::yes);

        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        TS_ASSERT( HeartConfig::Instance()->GetVisualizeWithVtk() );
        TS_ASSERT_EQUALS(HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().get().vtk(),
                         cp::yesno_type::yes);

        // Setting one doesn't change the others...
        TS_ASSERT( HeartConfig::Instance()->GetVisualizeWithCmgui() );
        TS_ASSERT( ! HeartConfig::Instance()->GetVisualizeWithMeshalyzer() );

        // Parameters file which does specify OutputVisualizer
        HeartConfig::Reset();
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");

        // Now the element exists...
        TS_ASSERT(HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().get().meshalyzer(),
                         cp::yesno_type::no);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().get().cmgui(),
                         cp::yesno_type::no);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->mpUserParameters->Simulation().get().OutputVisualizer().get().vtk(),
                         cp::yesno_type::yes);

        TS_ASSERT( ! HeartConfig::Instance()->GetVisualizeWithMeshalyzer() );
        TS_ASSERT( ! HeartConfig::Instance()->GetVisualizeWithCmgui() );
        TS_ASSERT( HeartConfig::Instance()->GetVisualizeWithVtk() );
    }

    void TestAdaptivityVariables()
    {
        // Defaults file doesn't have the AdaptivityParameters element
        TS_ASSERT( ! HeartConfig::Instance()->mpDefaultParameters->Numerical().AdaptivityParameters().present() );

        // Parameters file which doesn't specify AdaptivityParameters
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteEmpty.xml");

        // Read the user parameters directly - again element missing
        TS_ASSERT( ! HeartConfig::Instance()->mpUserParameters->Numerical().AdaptivityParameters().present() );

        // Get methods throw an exception
        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->GetTargetErrorForAdaptivity() );
        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->GetSigmaForAdaptivity() );
        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->GetMaxEdgeLengthForAdaptivity() );
        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->GetMinEdgeLengthForAdaptivity() );
        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->GetGradationForAdaptivity() );
        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->GetMaxNodesForAdaptivity() );
        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->GetNumberOfAdaptiveSweeps() );

        // Parameters file that does specify AdaptivityParameters
        HeartConfig::Reset();
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersFullFormat.xml");

        // Verify that the element is present
        TS_ASSERT( HeartConfig::Instance()->mpUserParameters->Numerical().AdaptivityParameters().present() );

        // Get methods
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetTargetErrorForAdaptivity(), 2.0, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetSigmaForAdaptivity(), 0.01, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMaxEdgeLengthForAdaptivity(), 0.04, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMinEdgeLengthForAdaptivity(), 0.005, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetGradationForAdaptivity(), 1.3, 1e-6 );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetMaxNodesForAdaptivity(), 1000u );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetNumberOfAdaptiveSweeps(), 5u );

        // Set methods
        HeartConfig::Instance()->SetAdaptivityParameters(0.1, 0.5, 0.3, 0.01, 2.0, 100, 7);

        TS_ASSERT_DELTA( HeartConfig::Instance()->GetTargetErrorForAdaptivity(), 0.1, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetSigmaForAdaptivity(), 0.5, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMaxEdgeLengthForAdaptivity(), 0.3, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMinEdgeLengthForAdaptivity(), 0.01, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetGradationForAdaptivity(), 2.0, 1e-6 );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetMaxNodesForAdaptivity(), 100u );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetNumberOfAdaptiveSweeps(), 7u );

        HeartConfig::Instance()->SetTargetErrorForAdaptivity(1.0);
        HeartConfig::Instance()->SetSigmaForAdaptivity(0.02);
        HeartConfig::Instance()->SetMaxEdgeLengthForAdaptivity(0.1);
        HeartConfig::Instance()->SetMinEdgeLengthForAdaptivity(0.001);
        HeartConfig::Instance()->SetGradationForAdaptivity(1.5);
        HeartConfig::Instance()->SetMaxNodesForAdaptivity(1000000u);
        HeartConfig::Instance()->SetNumberOfAdaptiveSweeps(3);

        TS_ASSERT_DELTA( HeartConfig::Instance()->GetTargetErrorForAdaptivity(), 1.0, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetSigmaForAdaptivity(), 0.02, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMaxEdgeLengthForAdaptivity(), 0.1, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMinEdgeLengthForAdaptivity(), 0.001, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetGradationForAdaptivity(), 1.5, 1e-6 );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetMaxNodesForAdaptivity(), (unsigned)1e6 );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetNumberOfAdaptiveSweeps(), 3u );

        TS_ASSERT_THROWS_ANYTHING( HeartConfig::Instance()->SetMinEdgeLengthForAdaptivity(1.0); );

        // Verify that if AdaptivityParameters element is not present then we can create it
        HeartConfig::Reset();
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteEmpty.xml");

        TS_ASSERT( ! HeartConfig::Instance()->mpUserParameters->Numerical().AdaptivityParameters().present() );
        HeartConfig::Instance()->SetAdaptivityParameters(0.1, 0.5, 0.3, 0.01, 2.0, 100, 7);
        TS_ASSERT( HeartConfig::Instance()->mpUserParameters->Numerical().AdaptivityParameters().present() );

        TS_ASSERT_DELTA( HeartConfig::Instance()->GetTargetErrorForAdaptivity(), 0.1, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetSigmaForAdaptivity(), 0.5, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMaxEdgeLengthForAdaptivity(), 0.3, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetMinEdgeLengthForAdaptivity(), 0.01, 1e-6 );
        TS_ASSERT_DELTA( HeartConfig::Instance()->GetGradationForAdaptivity(), 2.0, 1e-6 );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetMaxNodesForAdaptivity(), 100u );
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetNumberOfAdaptiveSweeps(), 7u );
    }

};

#endif /*TESTHEARTCONFIG_HPP_*/
