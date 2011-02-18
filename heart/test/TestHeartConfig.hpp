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
#include "ChasteEllipsoid.hpp"
#include "ChastePoint.hpp"
#include "Version.hpp"
#include "TetrahedralMesh.hpp"
#include "HeartFileFinder.hpp"
#include "Warnings.hpp"

class TestHeartConfig : public CxxTest::TestSuite
{
private:
    void setUp()
    {
        HeartConfig::Reset();
    }
    
public:
    void TestHeartConfigBasic() throw (Exception)
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

    void TestUserProvidedDifferentFromDefault() throw (Exception)
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

    void TestGetFunctions() throw (Exception)
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

        std::vector<boost::shared_ptr<AbstractStimulusFunction> > stimuli_applied;
        std::vector<ChasteCuboid<3> > stimulated_areas;
        HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_areas);

        TS_ASSERT_EQUALS(stimuli_applied.size(), 4u);
        TS_ASSERT_EQUALS(stimulated_areas.size(), 4u);

        TS_ASSERT_EQUALS(stimuli_applied[0]->GetStimulus(0), -25500.0);
        TS_ASSERT_EQUALS(stimuli_applied[0]->GetStimulus(0.6), 0.0);

        TS_ASSERT_EQUALS(stimuli_applied[2]->GetStimulus(0.0), 0.0);
        TS_ASSERT_EQUALS(stimuli_applied[2]->GetStimulus(2.0), -25500.0);
        TS_ASSERT_EQUALS(stimuli_applied[2]->GetStimulus(3.6), 0.0);

        //covering the 2D case
        std::vector<ChasteCuboid<2> > stimulated_areas_2D;
        std::vector<boost::shared_ptr<AbstractStimulusFunction> > stimuli_applied_2D;
        HeartConfig::Instance()->GetStimuli(stimuli_applied_2D, stimulated_areas_2D);

        TS_ASSERT_EQUALS(stimuli_applied_2D.size(), 4u);
        TS_ASSERT_EQUALS(stimulated_areas_2D.size(), 4u);

        //covering the 1D case
        std::vector<ChasteCuboid<1> > stimulated_areas_1D;
        std::vector<boost::shared_ptr<AbstractStimulusFunction> > stimuli_applied_1D;
        HeartConfig::Instance()->GetStimuli(stimuli_applied_1D, stimulated_areas_1D);

        TS_ASSERT_EQUALS(stimuli_applied_1D.size(), 4u);
        TS_ASSERT_EQUALS(stimulated_areas_1D.size(), 4u);

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
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "bjacobi")==0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);

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
        
        TS_ASSERT(HeartConfig::Instance()->IsPseudoEcgCalculationRequested());
        std::vector<ChastePoint<3> > pseudo_ecg_parameters;
        HeartConfig::Instance()->GetPseudoEcgElectrodePositions(pseudo_ecg_parameters);
        TS_ASSERT_EQUALS(pseudo_ecg_parameters.size(), 2u);
        for (unsigned dim=0; dim<3u; dim++)
        {
            TS_ASSERT_EQUALS(pseudo_ecg_parameters[0][dim], (double)dim);
            TS_ASSERT_EQUALS(pseudo_ecg_parameters[1][dim], (double)dim - 10.0);
        }

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), true);

        
        bool ground_second_electrode;
        unsigned axis_index;
        double magnitude, start_time, duration;
        HeartConfig::Instance()->GetElectrodeParameters(ground_second_electrode, axis_index, magnitude, start_time, duration);
        TS_ASSERT_EQUALS(ground_second_electrode, true);    
        TS_ASSERT_EQUALS(axis_index, 2u);    
        TS_ASSERT_EQUALS(magnitude, -11000.0);    
        TS_ASSERT_EQUALS(start_time, 1.0);    
        TS_ASSERT_EQUALS(duration, 2.0);
       
        
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

    void TestGetHeterogeneities() throw (Exception)
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
        std::vector<std::map<std::string, double> > parameter_settings;

        TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularTransmuralHeterogeneitiesRequested(), false);

        HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas_3D,
                                                        scale_factor_gks_3D,
                                                        scale_factor_ito_3D,
                                                        scale_factor_gkr_3D,
                                                        &parameter_settings);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularHeterogeneitiesSpecifiedByCuboids(), true);
        TS_ASSERT(cell_heterogeneity_areas_3D[0]->DoesContain(ChastePoint<3>(-1.0, 0, 0)));
        TS_ASSERT_EQUALS(scale_factor_gks_3D[0], 0.462);
        TS_ASSERT_EQUALS(scale_factor_ito_3D[0], 0.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_3D[0], 1.0);
        // The old scale factor elements have been removed from region 1, so values default to 1.0
        TS_ASSERT_EQUALS(scale_factor_gks_3D[1], 1.0);
        TS_ASSERT_EQUALS(scale_factor_ito_3D[1], 1.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_3D[1], 1.0);
        
        TS_ASSERT_EQUALS(parameter_settings[0].size(), 5u);
        TS_ASSERT_EQUALS(parameter_settings[1].size(), 3u);
        // Parameters get returned by name order since they're in a map
        std::map<std::string, double>::iterator param_it = parameter_settings[0].begin();
        TS_ASSERT_EQUALS(param_it->first, "ScaleFactorGkr");
        TS_ASSERT_EQUALS(param_it->second, 1.0);
        param_it++;
        TS_ASSERT_EQUALS(param_it->first, "ScaleFactorGks");
        TS_ASSERT_EQUALS(param_it->second, 0.462);
        param_it++;
        TS_ASSERT_EQUALS(param_it->first, "ScaleFactorIto");
        TS_ASSERT_EQUALS(param_it->second, 0.0);
        param_it++;
        TS_ASSERT_EQUALS(param_it->first, "example");
        TS_ASSERT_EQUALS(param_it->second, 0.0);
        param_it++;
        TS_ASSERT_EQUALS(param_it->first, "example2");
        TS_ASSERT_EQUALS(param_it->second, 2.0);

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
                                                        scale_factor_gkr_2D,
                                                        NULL);

        TS_ASSERT(cell_heterogeneity_areas_2D[0]->DoesContain(ChastePoint<2>(-1.0, 0)));
        TS_ASSERT_EQUALS(scale_factor_gks_2D[0], 0.462);
        TS_ASSERT_EQUALS(scale_factor_ito_2D[0], 0.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_2D[0], 1.0);
        TS_ASSERT_EQUALS(scale_factor_gks_2D[1], 1.0);
        TS_ASSERT_EQUALS(scale_factor_ito_2D[1], 1.0);
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
                                                        scale_factor_gkr_1D,
                                                        NULL);

        TS_ASSERT(cell_heterogeneity_areas_1D[0]->DoesContain(ChastePoint<1>(-1.0)));
        TS_ASSERT_EQUALS(scale_factor_gks_1D[0], 0.462);
        TS_ASSERT_EQUALS(scale_factor_ito_1D[0], 0.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_1D[0], 1.0);
        TS_ASSERT_EQUALS(scale_factor_gks_1D[1], 1.0);
        TS_ASSERT_EQUALS(scale_factor_ito_1D[1], 1.0);
        TS_ASSERT_EQUALS(scale_factor_gkr_1D[1], 1.0);

        for (unsigned i = 0; i<cell_heterogeneity_areas_1D.size();i++)
        {
            delete cell_heterogeneity_areas_1D[i];
        }

        ///////////////
        //Conductivity heterogeneities
        //////////////

        // Cuboid

        std::vector<AbstractChasteRegion<3>* > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        TS_ASSERT(conductivities_heterogeneity_areas[0]->DoesContain(ChastePoint<3>(1.95, 0, 0)));
        TS_ASSERT_EQUALS(intra_h_conductivities[0][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities[0][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities[1][0], 0.75);
        
        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
        {
            delete conductivities_heterogeneity_areas[region_index];
        }
        
        //cover the 2D case
        std::vector<AbstractChasteRegion<2>* > conductivities_heterogeneity_areas_2D;
        std::vector< c_vector<double,3> > intra_h_conductivities_2D;
        std::vector< c_vector<double,3> > extra_h_conductivities_2D;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_2D,
                                                                intra_h_conductivities_2D,
                                                                extra_h_conductivities_2D);
        
        TS_ASSERT(conductivities_heterogeneity_areas_2D[0]->DoesContain(ChastePoint<2>(1.95, 0)));
        TS_ASSERT_EQUALS(intra_h_conductivities_2D[0][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities_2D[0][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities_2D[1][0], 0.75);

        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas_2D.size(); region_index++)
        {
            delete conductivities_heterogeneity_areas_2D[region_index];
        }

        //cover the 1D case
        std::vector<AbstractChasteRegion<1>* > conductivities_heterogeneity_areas_1D;
        std::vector< c_vector<double,3> > intra_h_conductivities_1D;
        std::vector< c_vector<double,3> > extra_h_conductivities_1D;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_1D,
                                                                intra_h_conductivities_1D,
                                                                extra_h_conductivities_1D);
        
        TS_ASSERT(conductivities_heterogeneity_areas_1D[0]->DoesContain(ChastePoint<1>(1.95)));
        TS_ASSERT_EQUALS(intra_h_conductivities_1D[0][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities_1D[0][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities_1D[1][0], 0.75);

        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas_1D.size(); region_index++)
        {
            delete conductivities_heterogeneity_areas_1D[region_index];
        }

        // Ellipsoid

        std::vector<AbstractChasteRegion<3>* > conductivities_heterogeneity_areas_ellipsoid;
        std::vector< c_vector<double,3> > intra_h_conductivities_ellipsoid;
        std::vector< c_vector<double,3> > extra_h_conductivities_ellipsoid;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_ellipsoid,
        		intra_h_conductivities_ellipsoid,
        		extra_h_conductivities_ellipsoid);

        TS_ASSERT(conductivities_heterogeneity_areas_ellipsoid[3]->DoesContain(ChastePoint<3>(1, 0, 0)));
        TS_ASSERT_EQUALS(intra_h_conductivities_ellipsoid[3][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities_ellipsoid[3][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities_ellipsoid[4][0], 0.75);

        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas_ellipsoid.size(); region_index++)
        {
        	delete conductivities_heterogeneity_areas_ellipsoid[region_index];
        }

        //cover the 2D case
        std::vector<AbstractChasteRegion<2>* > conductivities_heterogeneity_areas_2D_ellipsoid;
        std::vector< c_vector<double,3> > intra_h_conductivities_2D_ellipsoid;
        std::vector< c_vector<double,3> > extra_h_conductivities_2D_ellipsoid;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_2D_ellipsoid,
        		intra_h_conductivities_2D_ellipsoid,
        		extra_h_conductivities_2D_ellipsoid);

        TS_ASSERT(conductivities_heterogeneity_areas_2D_ellipsoid[3]->DoesContain(ChastePoint<2>(1, 0)));
        TS_ASSERT_EQUALS(intra_h_conductivities_2D_ellipsoid[3][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities_2D_ellipsoid[3][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities_2D_ellipsoid[4][0], 0.75);

        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas_2D_ellipsoid.size(); region_index++)
        {
        	delete conductivities_heterogeneity_areas_2D_ellipsoid[region_index];
        }

        //cover the 1D case
        std::vector<AbstractChasteRegion<1>* > conductivities_heterogeneity_areas_1D_ellipsoid;
        std::vector< c_vector<double,3> > intra_h_conductivities_1D_ellipsoid;
        std::vector< c_vector<double,3> > extra_h_conductivities_1D_ellipsoid;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_1D_ellipsoid,
        		intra_h_conductivities_1D_ellipsoid,
        		extra_h_conductivities_1D_ellipsoid);

        TS_ASSERT(conductivities_heterogeneity_areas_1D_ellipsoid[3]->DoesContain(ChastePoint<1>(1)));
        TS_ASSERT_EQUALS(intra_h_conductivities_1D_ellipsoid[3][0], 2.75);
        TS_ASSERT_EQUALS(extra_h_conductivities_1D_ellipsoid[3][0], 8.0);
        TS_ASSERT_EQUALS(intra_h_conductivities_1D_ellipsoid[4][0], 0.75);

        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas_1D_ellipsoid.size(); region_index++)
        {
        	delete conductivities_heterogeneity_areas_1D_ellipsoid[region_index];
        }

        //extracellular conductivities
        c_vector<double, 3> extra_conductivities;
        HeartConfig::Instance()->GetExtracellularConductivities(extra_conductivities);
        TS_ASSERT_EQUALS(extra_h_conductivities[1][0], extra_conductivities[0]);

        TS_ASSERT_EQUALS(extra_conductivities[0], 7.0);
        TS_ASSERT_EQUALS(extra_conductivities[1], 7.0);
        TS_ASSERT_EQUALS(extra_conductivities[2], 7.0);
    }

    void TestIsMeshProvided() throw (Exception)
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

    void TestTransmuralHeterogeneities() throw (Exception)
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
                                                            scale_factor_gkr,
                                                            NULL);

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
            std::vector<AbstractChasteRegion<3>* > conductivities_heterogeneity_areas;
            std::vector< c_vector<double,3> > intra_h_conductivities;
            std::vector< c_vector<double,3> > extra_h_conductivities;
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                                          intra_h_conductivities,
                                                                                          extra_h_conductivities),
                                  "Definition of transmural layers is not allowed for conductivities heterogeneities, you may use fibre orientation support instead");

             //covering the case when the user specify transmural layers for stimulated areas (not yet supported)...
            std::vector<boost::shared_ptr<AbstractStimulusFunction> > stimuli_applied;
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
                                                                                  scale_factor_gkr,
                                                                                  NULL),
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
                                                                                  scale_factor_gkr,
                                                                                  NULL),
                                  "Fractions must be positive");

        }
        //covers the case when the user supplies negative numbers
        //this second test of supplied negative values is for coverage. As we first check that the summation is 1, it was impossible to have 3 negative values in one file.
        {
            HeartConfig::Reset();
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersCellHeterogeneities_negative_2.xml");

            std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
            std::vector<double> scale_factor_gks;
            std::vector<double> scale_factor_ito;
            std::vector<double> scale_factor_gkr;

            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                                                  scale_factor_gks,
                                                                                  scale_factor_ito,
                                                                                  scale_factor_gkr,
                                                                                  NULL),
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
                                                                                  scale_factor_gkr,
                                                                                  NULL),
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

            TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularHeterogeneitiesSpecifiedByCuboids(), false);
            TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                                                  scale_factor_gks,
                                                                                  scale_factor_ito,
                                                                                  scale_factor_gkr,
                                                                                  NULL),
                                  "Specification of cellular heterogeneities by cuboids and layers at the same time is not yet supported");

            TS_ASSERT_EQUALS(HeartConfig::Instance()->AreCellularHeterogeneitiesSpecifiedByCuboids(), true);
        }
    }
    
    void Test2dProblems() throw (Exception)
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

    void Test1dProblems() throw (Exception)
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


    void TestSetFunctions() throw (Exception)
    {
        // Start with a file that doesn't have much in it
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteEmpty.xml");
        
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

        std::vector<AbstractChasteRegion<3>* > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas.size(), 2u);
        TS_ASSERT_EQUALS(intra_h_conductivities.size(), 2u);
        TS_ASSERT_EQUALS(extra_h_conductivities.size(), 2u);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas[0]->DoesContain(ChastePoint<3>(0.0, 0.0, 0.0)), true);
        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas[0]->DoesContain(ChastePoint<3>(-1.5, -1.5, -1.5)), false);
        TS_ASSERT_EQUALS(intra_h_conductivities[0][0], 2.5);
        TS_ASSERT_EQUALS(extra_h_conductivities[0][0], 8.5);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas[1]->DoesContain(ChastePoint<3>(-1.5, -1.5, -1.5)), true);
        TS_ASSERT_EQUALS(intra_h_conductivities[1][0], 1.0);
        
        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
        {
            delete conductivities_heterogeneity_areas[region_index];
        }
        
        // Test ellipsoid
        std::vector< c_vector<double,3> > intra_conductivities_ellipsoid;
        std::vector< c_vector<double,3> > extra_conductivities_ellipsoid;

        std::vector<ChasteEllipsoid<3> > input_areas_ellipsoid;
        ChastePoint<3> centre_1(-1.0, -1.0, -1.0);
        ChastePoint<3> radii_1( 1.5,  2.0,  3.0);
        input_areas_ellipsoid.push_back(ChasteEllipsoid<3> (centre_1, radii_1));
        intra_conductivities_ellipsoid.push_back( Create_c_vector(2.5, 2.5, 2.5) );
        extra_conductivities_ellipsoid.push_back( Create_c_vector(8.5, 8.5, 8.5) );

        ChastePoint<3> centre_2(-2.0, -2.0, -2.0);
        ChastePoint<3> radii_2(1.0,  2.0,  3.0);
        input_areas_ellipsoid.push_back(ChasteEllipsoid<3> (centre_2, radii_2));
        intra_conductivities_ellipsoid.push_back( Create_c_vector(1.0, 0.5, 0.4) );
        extra_conductivities_ellipsoid.push_back( Create_c_vector(7.0, 6.5, 6.4) );

        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(input_areas_ellipsoid, intra_conductivities_ellipsoid, extra_conductivities_ellipsoid);

        TS_ASSERT(HeartConfig::Instance()->GetConductivityHeterogeneitiesProvided());

        std::vector<AbstractChasteRegion<3>* > conductivities_heterogeneity_areas_ellipsoid;
        std::vector< c_vector<double,3> > intra_h_conductivities_ellipsoid;
        std::vector< c_vector<double,3> > extra_h_conductivities_ellipsoid;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas_ellipsoid,
                                                                intra_h_conductivities_ellipsoid,
                                                                extra_h_conductivities_ellipsoid);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas_ellipsoid.size(), 2u);
        TS_ASSERT_EQUALS(intra_h_conductivities_ellipsoid.size(), 2u);
        TS_ASSERT_EQUALS(extra_h_conductivities_ellipsoid.size(), 2u);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas_ellipsoid[0]->DoesContain(ChastePoint<3>(0.0, 0.0, 0.0)), true);
        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas_ellipsoid[0]->DoesContain(ChastePoint<3>(-2.5, -1.5, -1.5)), false);
        TS_ASSERT_EQUALS(intra_h_conductivities_ellipsoid[0][0], 2.5);
        TS_ASSERT_EQUALS(extra_h_conductivities_ellipsoid[0][0], 8.5);

        TS_ASSERT_EQUALS(conductivities_heterogeneity_areas_ellipsoid[1]->DoesContain(ChastePoint<3>(-1.5, -1.5, -1.5)), true);
        TS_ASSERT_EQUALS(intra_h_conductivities_ellipsoid[1][0], 1.0);
       
        // freeing memory allcated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas_ellipsoid.size(); region_index++)
        {
            delete conductivities_heterogeneity_areas_ellipsoid[region_index];
        }
        
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

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), false);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), true);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), false);
        
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

        std::map<unsigned, double> multiple_bath_conductivities;        
        multiple_bath_conductivities[2] = 3.0;
        multiple_bath_conductivities[4] = 4.0;
        multiple_bath_conductivities[5] = 5.0;

        HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetBathConductivity(), 150);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetBathConductivity(2), 3.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetBathConductivity(4), 4.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetBathConductivity(5), 5.0);

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetBathConductivity(0), "Region label 0 is reserved for tissue");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetBathConductivity(3), "Bath conductivity not defined for region 3");

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
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdeTimeStep(0.2), "Ode time-step should not be greater than PDE time-step");

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(-0.1, 0.1, 0.1), "Ode time-step should be positive");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, -0.1, 0.1), "Pde time-step should be positive");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, -0.1), "Printing time-step should be positive");

        //Throws when we try to print more often than the pde time step
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.1), "Printing time-step should not be smaller than PDE time-step");
         //Throws when printing step is not a multiple of pde time step
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.3), "Printing time-step should be a multiple of PDE time step");

        //Shouldn't throw because 10 is a multiple of 0.1
        TS_ASSERT_THROWS_NOTHING(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, 10.0));

        // Throws because ode time step is bigger than pde time step
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.3), "Printing time-step should be a multiple of PDE time step");

        /*
         *  Set ODE, PDE, and printing timestep to something meaningful and test SetCheckpointSimulation() exceptions.
         */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.2, 0.4);
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetCheckpointSimulation(true, 1.0, 3), "Checkpoint time-step should be a multiple of printing time-step");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetCheckpointSimulation(true, -2.0, 3), "Checkpoint time-step should be positive");

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

        HeartConfig::Instance()->SetKSPSolver("chebychev");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "chebychev")==0);

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetKSPSolver("foobar"),"Unknown solver type provided");

        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "jacobi")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "bjacobi")==0);

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

        HeartConfig::Instance()->SetKSPPreconditioner("twolevelsblockdiagonal");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "twolevelsblockdiagonal")==0);

        HeartConfig::Instance()->SetKSPPreconditioner("none");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "none")==0);

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetKSPPreconditioner("foobar"),
                "Unknown preconditioner type provided");

        // Mesh partitioning method
        HeartConfig::Instance()->SetMeshPartitioning("dumb");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::DUMB);
        HeartConfig::Instance()->SetMeshPartitioning("metis");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
        HeartConfig::Instance()->SetMeshPartitioning("parmetis");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
        HeartConfig::Instance()->SetMeshPartitioning("petsc");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION);

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetMeshPartitioning("magic"),
                              "Unknown mesh partitioning method provided");

        // Tests for set functions of postprocessing
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingSectionPresent(), true); // It's in the defaults, but empty
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingRequested(), false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsApdMapsRequested(), false);
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
        
        TS_ASSERT(!HeartConfig::Instance()->IsPseudoEcgCalculationRequested());
        std::vector<ChastePoint<3> > pseudo_ecg_parameters, pseudo_ecg_parameters_get;
        ChastePoint<3> electrode_point(0.0, 1.5, -2.5);
        pseudo_ecg_parameters.push_back(electrode_point);
        HeartConfig::Instance()->SetPseudoEcgElectrodePositions(pseudo_ecg_parameters);
        TS_ASSERT(HeartConfig::Instance()->IsPseudoEcgCalculationRequested());
        HeartConfig::Instance()->GetPseudoEcgElectrodePositions(pseudo_ecg_parameters_get);
        TS_ASSERT_EQUALS(pseudo_ecg_parameters_get.size(), 1u);
        for (unsigned dim=0; dim<3u; dim++)
        {
            TS_ASSERT_EQUALS(pseudo_ecg_parameters_get[0][dim], electrode_point[dim]);
        }
        
        
        
        bool ground_second_electrode;
        unsigned axis_index;
        double magnitude, start_time, duration;
        //The electrodes section does not exist yet.
        TS_ASSERT_THROWS_THIS( HeartConfig::Instance()->GetElectrodeParameters(
                               ground_second_electrode, axis_index, magnitude, start_time, duration),
                               "Attempted to get electrodes that have not been defined.");
        HeartConfig::Instance()->SetElectrodeParameters(false, 2, 1066, 0.5, 0.5);
        HeartConfig::Instance()->GetElectrodeParameters(ground_second_electrode, axis_index, magnitude, start_time, duration);
        TS_ASSERT_EQUALS(ground_second_electrode, false);    
        TS_ASSERT_EQUALS(axis_index, 2u);    
        TS_ASSERT_EQUALS(magnitude, 1066.0);    
        TS_ASSERT_EQUALS(start_time, 0.5);    
        TS_ASSERT_EQUALS(duration, 0.5);
        
        // This is a temporary internal boolean until we're happy that users can be let loose on the functionality!
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseStateVariableInterpolation(), false);
        HeartConfig::Instance()->SetUseStateVariableInterpolation();
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseStateVariableInterpolation(), true);
        HeartConfig::Instance()->SetUseStateVariableInterpolation(false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseStateVariableInterpolation(), false);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseMassLumping(), false);
        HeartConfig::Instance()->SetUseMassLumping();
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseMassLumping(), true);
        HeartConfig::Instance()->SetUseMassLumping(false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseMassLumping(), false);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseReactionDiffusionOperatorSplitting(), false);
        HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplitting();
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseReactionDiffusionOperatorSplitting(), true);
        HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplitting(false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseReactionDiffusionOperatorSplitting(), false);

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseMassLumpingForPrecond(), false);
        HeartConfig::Instance()->SetUseMassLumpingForPrecond();
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseMassLumpingForPrecond(), true);
        HeartConfig::Instance()->SetUseMassLumpingForPrecond(false);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetUseMassLumpingForPrecond(), false);

    }

    void TestWrite() throw (Exception)
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
        
        // If we Write when the schema fixed location doesn't exist, and the schema isn't in
        // the current directory, we should get a warning.  This could occur if an executable
        // user doesn't have the corresponding source tree.
        HeartConfig::Reset();
        FileFinder::FakePath(RelativeTo::ChasteSourceRoot, "/not-there");
        try
        {
            output_file_handler.SetArchiveDirectory();
            HeartConfig::Instance()->Write(true);
            FileFinder::StopFaking();
        }
        catch (...)
        {
            FileFinder::StopFaking();
            throw;
        }
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),
                         "Unable to locate schema file ChasteParameters_2_2.xsd. You will need to ensure it is available when resuming from the checkpoint.");
    }

    void TestArchiving() throw (Exception)
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
            HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersResumeSimulationFullFormat.xml");
            input_arch >> (*p_heart_config);
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 20.0);
            TS_ASSERT(p_heart_config->GetDefaultIonicModel().Hardcoded().present());
            TS_ASSERT_EQUALS( user_ionic, p_heart_config->GetDefaultIonicModel().Hardcoded().get());

            // Check that the resume parameters have overridden everything they should have
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetOdeTimeStep(), 0.01, 1e-12);
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPdeTimeStep(), 0.01, 1e-12);
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), 0.01, 1e-12);
            TS_ASSERT(HeartConfig::Instance()->GetUseAbsoluteTolerance());
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetAbsoluteTolerance(), 1e-5, 1e-12);
            TS_ASSERT_EQUALS(strcmp(HeartConfig::Instance()->GetKSPSolver(), "cg"), 0);
            TS_ASSERT_EQUALS(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "none"), 0);
            TS_ASSERT(HeartConfig::Instance()->IsAdaptivityParametersPresent());
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetTargetErrorForAdaptivity(), 0.0);
            TS_ASSERT(HeartConfig::Instance()->IsPostProcessingRequested());
            TS_ASSERT(HeartConfig::Instance()->IsApdMapsRequested());
            std::vector<std::pair<double,double> > apd_maps;
            HeartConfig::Instance()->GetApdMaps(apd_maps);
            TS_ASSERT_EQUALS(apd_maps.size(), 1u);
            TS_ASSERT_DELTA(apd_maps[0].first, 70.0, 1e-12);
            TS_ASSERT_DELTA(apd_maps[0].second, -20.0, 1e-12);
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
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteInconsistent.xml"),
                "Ode time-step should not be greater than PDE time-step");
        HeartConfig::Reset();
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteInconsistentCheckpointTimestep.xml"),
                "Checkpoint time-step should be a multiple of printing time-step");

        // Stimulus defined with end time but not period
        HeartConfig::Reset();
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersWrongSimpleStimulusDefinition.xml");

        std::vector<boost::shared_ptr<AbstractStimulusFunction> > stimuli_applied;
        std::vector<ChasteCuboid<3> > stimulated_areas;        

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_areas),
                "Stop time can not be defined for SimpleStimulus. Use Duration instead.");

        //Can't open a directory
        HeartConfig::Reset();
        HeartConfig::Instance()->SetOutputDirectory("../../../");
        TS_ASSERT_THROWS_CONTAINS(HeartConfig::Instance()->Write(), "due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");

        // Can't open a file for writing
        std::string command = OutputFileHandler::GetChasteTestOutputDirectory() + "no_write_access";
        mkdir(command.c_str(), 0);
        chmod(command.c_str(), 0); // in case the directory already exists
        HeartConfig::Instance()->SetOutputDirectory("no_write_access");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->Write(false, ""),
                              "Could not open XML file in HeartConfig");
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
    void TestVersioning() throw (Exception)
    {
        // Test we can recognise known versions
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetVersionFromNamespace(""), 1001u);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetVersionFromNamespace("https://chaste.comlab.ox.ac.uk/nss/parameters/2_0"), 2000u);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetVersionFromNamespace("https://chaste.comlab.ox.ac.uk/nss/parameters/2_1"), 2001u);
        // and exceptions
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetVersionFromNamespace("https://chaste.comlab.ox.ac.uk/nss/parameters/1__1"),
                              "https://chaste.comlab.ox.ac.uk/nss/parameters/1__1 is not a recognised Chaste parameters namespace.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetVersionFromNamespace("bob"),
                              "bob is not a recognised Chaste parameters namespace.");
        
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
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingSectionPresent(), false); // Comes from latest defaults
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

        // Check that release 2.0 xml can be loaded with release 2.0 schema
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease2_0.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease2_0.xml");

        // Check that release 2.0 xml can be loaded with latest schema
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease2_0.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease2_0.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingSectionPresent(), true);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->IsPostProcessingRequested(), false);
        
        // We removed ilu in version 2.1; throw an exception if an older parameters file uses it
        HeartConfig::Reset();
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ilu_preconditioner.xml"),
                              "PETSc does not have a parallel implementation of ilu, so we no longer allow it as an option.  Use bjacobi instead.");

        // Check that release 2.1 xml can be loaded with release 2.1 schema
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease2_1.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);

        // Check that release 2.1 xml can be loaded with latest schema
        HeartConfig::Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease2_1.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::METIS_LIBRARY);
    }

    /**
     * Test whether we can use a schema that has spaces in the pathname.
     * This gives some indication of whether Chaste will cope being checked out into
     * a path with spaces.
     */
    void TestSpacesInPath() throw (Exception)
    {
        HeartConfig::Reset();
        HeartConfig::SchemaLocationsMap schema_locations;
        schema_locations[""] = std::string(ChasteBuildInfo::GetRootDir()) + "/heart/test/data/xml/schema with spaces.xsd";
        HeartConfig::Instance()->SetFixedSchemaLocations(schema_locations);
        HeartConfig::Instance()->SetDefaultsFile("heart/test/data/xml/ChasteDefaultsRelease1_1.xml");
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersRelease1_1.xml");
    }

    void TestGetOuputVariablesFromXML() throw (Exception)
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

    void TestSetAndGetOuputVariables() throw (Exception)
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

    void TestSetAndGetArchivingStuff() throw (Exception)
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

        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/ChasteParametersResumeSimulationFullFormat.xml");
        TS_ASSERT(!HeartConfig::Instance()->IsSimulationDefined());
        TS_ASSERT(HeartConfig::Instance()->IsSimulationResumed());

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetArchivedSimulationDir().GetAbsolutePath(),
                         OutputFileHandler::GetChasteTestOutputDirectory() + "ChasteResults_10ms");
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

        std::vector<boost::shared_ptr<AbstractStimulusFunction> > stimuli_applied;
        std::vector<ChasteCuboid<3> > stimulated_area;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_area), "Stimuli information is not available in a resumed simulation.")

        std::vector<AbstractChasteRegion<3>* > cell_heterogeneity_areas;
        std::vector<double> scale_factor_gks;
        std::vector<double> scale_factor_ito;
        std::vector<double> scale_factor_gkr;
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas, scale_factor_gks, scale_factor_ito, scale_factor_gkr, NULL),
                              "CellHeterogeneities information is not available in a resumed simulation.");
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->GetConductivityHeterogeneitiesProvided(),
                              "ConductivityHeterogeneities information is not available in a resumed simulation.");

        for (unsigned i = 0; i<cell_heterogeneity_areas.size();i++)
        {
            delete cell_heterogeneity_areas[i];
        }

        std::vector<AbstractChasteRegion<3>* > conductivitiesHeterogeneityAreas;
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

    void TestOutputVisualizerSettings() throw (Exception)
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

    void TestAdaptivityVariables() throw (Exception)
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
