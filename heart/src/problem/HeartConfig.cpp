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


#include "HeartConfig.hpp"
using namespace xsd::cxx::tree;

std::auto_ptr<HeartConfig> HeartConfig::mpInstance;

HeartConfig* HeartConfig::Instance()
{
    if (mpInstance.get() == NULL)
    {
        mpInstance.reset(new HeartConfig);
    }
    return mpInstance.get();
}

HeartConfig::HeartConfig()
{
    assert(mpInstance.get() == NULL);
    mpDefaultParameters = NULL;
    mpUserParameters = NULL;
    mpDefaultParameters = ReadFile("ChasteDefaults.xml");
    mpUserParameters = mpDefaultParameters;
    //CheckTimeSteps(); // necessity of this line of code is not tested -- remove with caution!
}

HeartConfig::~HeartConfig()
{
    if (mpUserParameters != mpDefaultParameters)
    {
        delete mpUserParameters;
    }

    delete mpDefaultParameters;    
}

void HeartConfig::SetDefaultsFile(std::string fileName)
{
    bool same_target = (mpUserParameters == mpDefaultParameters);
    
    delete mpDefaultParameters;
    mpDefaultParameters = ReadFile(fileName);

    if (same_target)
    {
        mpUserParameters = mpDefaultParameters;
    }
    CheckTimeSteps();
}

void HeartConfig::Write(std::string dirName, std::string fileName)
{
    //Output file
    OutputFileHandler output_file_handler(dirName, false);
    out_stream p_file = output_file_handler.OpenOutputFile(fileName);

    //Schema map
    //Note - this location is relative to where we are storing the xml
    xml_schema::namespace_infomap map;
    char buf[10000];
    std::string absolute_path_to_xsd=getcwd(buf, 10000);
    absolute_path_to_xsd += "/heart/src/io/ChasteParameters.xsd";
    map[""].schema = absolute_path_to_xsd;
    
    ChasteParameters(*p_file, *mpUserParameters, map);   
}
chaste_parameters_type* HeartConfig::ReadFile(std::string fileName)
{
    // get the parameters using the method 'ChasteParameters(filename)',
    // which returns a std::auto_ptr. We don't want to use a std::auto_ptr because
    // it will delete memory when out of scope, or no longer point when it is copied,
    // so we reallocate memory using a normal pointer and copy the data to there
    try
    {
        std::auto_ptr<chaste_parameters_type> p_default(ChasteParameters(fileName));
        return new chaste_parameters_type(*p_default);
    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         //More clunky memory management
         mpUserParameters = NULL;
         mpDefaultParameters = NULL;
         EXCEPTION("XML parsing error in configuration file: " + fileName);
    }
}

void HeartConfig::SetParametersFile(std::string fileName)
{
    // handles multiple calls to the method in the same context
    if (mpUserParameters != mpDefaultParameters)
    {        
        delete mpUserParameters;         
    }
    mpUserParameters = ReadFile(fileName);
}

chaste_parameters_type* HeartConfig::UserParameters()
{
    return mpUserParameters;
}

chaste_parameters_type* HeartConfig::DefaultParameters()
{
    return mpDefaultParameters;
}

void HeartConfig::Reset()
{
    //Throw it away
    mpInstance.reset(0);
    //Make a new one
    mpInstance.reset(new HeartConfig);
}

template<class TYPE>
TYPE* HeartConfig::DecideLocation(TYPE* ptr1, TYPE* ptr2, const std::string& nameParameter) const
{
    if (ptr1->present())
    {
        return ptr1;
    }
    if (ptr2->present())
    {
        return ptr2;
    }
    EXCEPTION("No " + nameParameter + " provided (neither default nor user defined)");
}

double HeartConfig::GetSimulationDuration() const
{
    return DecideLocation( & mpUserParameters->Simulation().SimulationDuration(),
                           & mpDefaultParameters->Simulation().SimulationDuration(),
                           "SimulationDuration")->get();
}

domain_type HeartConfig::GetDomain() const
{
    return DecideLocation( & mpUserParameters->Simulation().Domain(),
                           & mpDefaultParameters->Simulation().Domain(),
                           "Domain")->get();
}

ionic_model_type HeartConfig::GetIonicModel() const
{
    return DecideLocation( & mpUserParameters->Simulation().IonicModel(),
                           & mpDefaultParameters->Simulation().IonicModel(),
                           "IonicModel")->get();
}

bool HeartConfig::GetIsMeshProvided() const
{
	try
	{
    	DecideLocation( & mpUserParameters->Simulation().Mesh(),
                    	& mpDefaultParameters->Simulation().Mesh(),
                    	"Mesh");                    	
    	return true;
	}
	catch (Exception& e)
	{			
		return false;
	}
}	

bool HeartConfig::GetCreateSlab() const
{
    return (DecideLocation( & mpUserParameters->Simulation().Mesh(),
                            & mpDefaultParameters->Simulation().Mesh(),
                            "Mesh")->get().Slab().present());
}

bool HeartConfig::GetLoadMesh() const
{
	return (DecideLocation( & mpUserParameters->Simulation().Mesh(),
                        	& mpDefaultParameters->Simulation().Mesh(),
                        	"Mesh")->get().LoadMesh().present());
}
 
void HeartConfig::GetSlabDimensions(c_vector<double, 3>& slabDimensions) const
{
    assert(GetCreateSlab());
    
    optional<slab_type, false> slab_dimensions = DecideLocation( & mpUserParameters->Simulation().Mesh(),
                                                                  & mpDefaultParameters->Simulation().Mesh(),
                                                                  "Slab")->get().Slab();
    
    slabDimensions[0] = slab_dimensions->SlabX();
    slabDimensions[1] = slab_dimensions->SlabY();
    slabDimensions[2] = slab_dimensions->SlabZ();
}

double HeartConfig::GetInterNodeSpace() const
{
    assert(GetCreateSlab());

    return DecideLocation( & mpUserParameters->Simulation().Mesh(),
                           & mpDefaultParameters->Simulation().Mesh(),
                           "Slab")->get().Slab()->InterNodeSpace();
}

std::string HeartConfig::GetMeshName() const
{
    assert(GetLoadMesh());
    
    return DecideLocation( & mpUserParameters->Simulation().Mesh(),
                           & mpDefaultParameters->Simulation().Mesh(),
                           "LoadMesh")->get().LoadMesh()->name();
}

media_type HeartConfig::GetConductivityMedia() const
{
    assert(GetLoadMesh());

    return DecideLocation( & mpUserParameters->Simulation().Mesh(),
                           & mpDefaultParameters->Simulation().Mesh(),
                           "LoadMesh")->get().LoadMesh()->conductivity_media();           
}

void HeartConfig::GetStimuli(std::vector<SimpleStimulus>& stimuliApplied, std::vector<ChasteCuboid>& stimulatedAreas) const
{

    simulation_type::Stimuli::_xsd_Stimuli_::Stimuli::Stimulus::container&
         stimuli = DecideLocation( & mpUserParameters->Simulation().Stimuli(),
                           & mpDefaultParameters->Simulation().Stimuli(),
                           "Stimuli")->get().Stimulus();
    for (simulation_type::Stimuli::_xsd_Stimuli_::Stimuli::Stimulus::iterator i = stimuli.begin();
         i != stimuli.end();
         ++i)
    {
        stimulus_type stimulus(*i);
        point_type point_a = stimulus.Location().CornerA();
        point_type point_b = stimulus.Location().CornerB();

        ChastePoint<3> chaste_point_a ( point_a.x(),
                                        point_a.y(),
                                        point_a.z());

        ChastePoint<3> chaste_point_b ( point_b.x(),
                                        point_b.y(),
                                        point_b.z());

        stimuliApplied.push_back( SimpleStimulus(stimulus.Strength(), stimulus.Duration(), stimulus.Delay() ) );
        stimulatedAreas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
    }
}

void HeartConfig::GetCellHeterogeneities(std::vector<ChasteCuboid>& cellHeterogeneityAreas,
    								     std::vector<double>& scaleFactorGks,
    							         std::vector<double>& scaleFactorIto) const
{
    simulation_type::CellHeterogeneities::_xsd_CellHeterogeneities_::CellHeterogeneities::CellHeterogeneity::container&
         cell_heterogeneity = DecideLocation( & mpUserParameters->Simulation().CellHeterogeneities(),
                                                 & mpDefaultParameters->Simulation().CellHeterogeneities(),
                                                 "CellHeterogeneities")->get().CellHeterogeneity();

    for (simulation_type::CellHeterogeneities::_xsd_CellHeterogeneities_::CellHeterogeneities::CellHeterogeneity::iterator i = cell_heterogeneity.begin();
         i != cell_heterogeneity.end();
         ++i)
    {
        cell_heterogeneity_type ht(*i);
        point_type point_a = ht.Location().CornerA();
        point_type point_b = ht.Location().CornerB();

        ChastePoint<3> chaste_point_a (point_a.x(),
                                       point_a.y(),
                                       point_a.z());

        ChastePoint<3> chaste_point_b (point_b.x(),
                                       point_b.y(),
                                       point_b.z());

        scaleFactorGks.push_back (ht.ScaleFactorGks());
        scaleFactorIto.push_back (ht.ScaleFactorIto());
        cellHeterogeneityAreas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
    }
}

void HeartConfig::GetConductivityHeterogeneities(std::vector<ChasteCuboid>& conductivitiesHeterogeneityAreas,
					  			 	             std::vector< c_vector<double,3> >& intraConductivities,
									             std::vector< c_vector<double,3> >& extraConductivities) const
{
    simulation_type::ConductivityHeterogeneities::_xsd_ConductivityHeterogeneities_::ConductivityHeterogeneities::ConductivityHeterogeneity::container&
         conductivity_heterogeneity = DecideLocation( & mpUserParameters->Simulation().ConductivityHeterogeneities(),
                                                         & mpDefaultParameters->Simulation().ConductivityHeterogeneities(),
                                                           "CellHeterogeneities")->get().ConductivityHeterogeneity();

    for (simulation_type::ConductivityHeterogeneities::_xsd_ConductivityHeterogeneities_::ConductivityHeterogeneities::ConductivityHeterogeneity::iterator i = conductivity_heterogeneity.begin();
         i != conductivity_heterogeneity.end();
         ++i)
    {
        conductivity_heterogeneity_type ht(*i);
        point_type point_a = ht.Location().CornerA();
        point_type point_b = ht.Location().CornerB();

        ChastePoint<3> chaste_point_a (point_a.x(),
                                       point_a.y(),
                                       point_a.z());

        ChastePoint<3> chaste_point_b (point_b.x(),
                                       point_b.y(),
                                       point_b.z());

        conductivitiesHeterogeneityAreas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );

        if (ht.IntracellularConductivities().present())
        {
            double intra_x = ht.IntracellularConductivities().get().longi();
            double intra_y = ht.IntracellularConductivities().get().trans();
            double intra_z = ht.IntracellularConductivities().get().normal();

            intraConductivities.push_back( Create_c_vector(intra_x, intra_y, intra_z) );
        }
        else
        {
            c_vector<double, 3> intra_conductivities;
            GetIntracellularConductivities(intra_conductivities);
            intraConductivities.push_back(intra_conductivities);
        }

        if (ht.ExtracellularConductivities().present())
        {
            double extra_x = ht.ExtracellularConductivities().get().longi();
            double extra_y = ht.ExtracellularConductivities().get().trans();
            double extra_z = ht.ExtracellularConductivities().get().normal();

            extraConductivities.push_back( Create_c_vector(extra_x, extra_y, extra_z) );
        }
        else
        {
            c_vector<double, 3> extra_conductivities;
            GetExtracellularConductivities(extra_conductivities);
            extraConductivities.push_back(extra_conductivities);
        }

    }
}

std::string HeartConfig::GetOutputDirectory() const
{
    return DecideLocation( & mpUserParameters->Simulation().OutputDirectory(),
                           & mpDefaultParameters->Simulation().OutputDirectory(),
                           "OutputDirectory")->get();
}

std::string HeartConfig::GetOutputFilenamePrefix() const
{
    return DecideLocation( & mpUserParameters->Simulation().OutputFilenamePrefix(),
                           & mpDefaultParameters->Simulation().OutputFilenamePrefix(),
                           "OutputFilenamePrefix")->get();	
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 3>& intraConductivities) const
{
    optional<conductivities_type, false>* intra_conductivities  = DecideLocation( & mpUserParameters->Physiological().IntracellularConductivities(),
                                                                                  & mpDefaultParameters->Physiological().IntracellularConductivities(),
                                                                                  "IntracellularConductivities");
    double intra_x_cond = intra_conductivities->get().longi();
    double intra_y_cond = intra_conductivities->get().trans();
    double intra_z_cond = intra_conductivities->get().normal();;

    assert(intra_y_cond != DBL_MAX); 
    assert(intra_z_cond != DBL_MAX); 

    intraConductivities[0] = intra_x_cond;
    intraConductivities[1] = intra_y_cond;
    intraConductivities[2] = intra_z_cond;
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 2>& intraConductivities) const
{
    optional<conductivities_type, false>* intra_conductivities  = DecideLocation( & mpUserParameters->Physiological().IntracellularConductivities(),
                                                                                  & mpDefaultParameters->Physiological().IntracellularConductivities(),
                                                                                  "IntracellularConductivities");
    double intra_x_cond = intra_conductivities->get().longi();
    double intra_y_cond = intra_conductivities->get().trans();

    assert(intra_y_cond != DBL_MAX);  

    intraConductivities[0] = intra_x_cond;
    intraConductivities[1] = intra_y_cond;
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 1>& intraConductivities) const
{
    optional<conductivities_type, false>* intra_conductivities  = DecideLocation( & mpUserParameters->Physiological().IntracellularConductivities(),
                                                                                  & mpDefaultParameters->Physiological().IntracellularConductivities(),
                                                                                  "IntracellularConductivities");
    double intra_x_cond = intra_conductivities->get().longi();

    intraConductivities[0] = intra_x_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 3>& extraConductivities) const
{
    optional<conductivities_type, false>* extra_conductivities  = DecideLocation( & mpUserParameters->Physiological().ExtracellularConductivities(),
                                                                                  & mpDefaultParameters->Physiological().ExtracellularConductivities(),
                                                                                  "ExtracellularConductivities");
    double extra_x_cond = extra_conductivities->get().longi();
    double extra_y_cond = extra_conductivities->get().trans();
    double extra_z_cond = extra_conductivities->get().normal();;

    assert(extra_y_cond != DBL_MAX); 
    assert(extra_z_cond != DBL_MAX); 

    extraConductivities[0] = extra_x_cond;
    extraConductivities[1] = extra_y_cond;
    extraConductivities[2] = extra_z_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 2>& extraConductivities) const
{
    optional<conductivities_type, false>* extra_conductivities  = DecideLocation( & mpUserParameters->Physiological().ExtracellularConductivities(),
                                                                                  & mpDefaultParameters->Physiological().ExtracellularConductivities(),
                                                                                  "ExtracellularConductivities");
    double extra_x_cond = extra_conductivities->get().longi();
    double extra_y_cond = extra_conductivities->get().trans();

    assert(extra_y_cond != DBL_MAX);  

    extraConductivities[0] = extra_x_cond;
    extraConductivities[1] = extra_y_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 1>& extraConductivities) const
{
    optional<conductivities_type, false>* extra_conductivities  = DecideLocation( & mpUserParameters->Physiological().ExtracellularConductivities(),
                                                                                  & mpDefaultParameters->Physiological().ExtracellularConductivities(),
                                                                                  "ExtracellularConductivities");
    double extra_x_cond = extra_conductivities->get().longi();

    extraConductivities[0] = extra_x_cond;
}

double HeartConfig::GetSurfaceAreaToVolumeRatio() const
{
    /*surface area to volume ratio: 1/cm*/
    return DecideLocation( & mpUserParameters->Physiological().SurfaceAreaToVolumeRatio(),
                           & mpDefaultParameters->Physiological().SurfaceAreaToVolumeRatio(),
                           "SurfaceAreaToVolumeRatio")->get();
}

double HeartConfig::GetCapacitance() const
{
    //         capacitance                 : uF/cm^2
    return DecideLocation( & mpUserParameters->Physiological().Capacitance(),
                           & mpDefaultParameters->Physiological().Capacitance(),
                           "Capacitance")->get();
}

double HeartConfig::GetOdeTimeStep() const
{
    return DecideLocation( & mpUserParameters->Numerical().TimeSteps(),
                           & mpDefaultParameters->Numerical().TimeSteps(),
                           "ode TimeStep")->get().ode();
}

double HeartConfig::GetPdeTimeStep() const
{
    return DecideLocation( & mpUserParameters->Numerical().TimeSteps(),
                           & mpDefaultParameters->Numerical().TimeSteps(),
                           "pde TimeStep")->get().pde();
}

double HeartConfig::GetPrintingTimeStep() const
{
    return DecideLocation( & mpUserParameters->Numerical().TimeSteps(),
                           & mpDefaultParameters->Numerical().TimeSteps(),
                           "printing TimeStep")->get().printing();
}

bool HeartConfig::GetUseAbsoluteTolerance() const
{
     /*
      * Note that it may be the case that absolute tolerance exists in the default
      * parameters file, but has been overridden in the user parameters
      */
     if (mpUserParameters->Numerical().KSPTolerances().get().KSPRelative().present() )
     {
        return false;
     }

     return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                                             & mpDefaultParameters->Numerical().KSPTolerances(),
                                             "KSPTolerances")->get().KSPAbsolute().present();
}

double HeartConfig::GetAbsoluteTolerance() const
{
	if (!GetUseAbsoluteTolerance())
    {
        EXCEPTION("Absolute tolerance is not set in Chaste parameters");
    }	
    return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                           & mpDefaultParameters->Numerical().KSPTolerances(),
                           "KSPTolerances")->get().KSPAbsolute().get();
}

bool HeartConfig::GetUseRelativeTolerance() const
{
    /*
      * Note that it may be the case that relative tolerance exists in the default
      * parameters file, but has been overridden in the user parameters
      */
     
     
     if (mpUserParameters->Numerical().KSPTolerances().get().KSPAbsolute().present() )
     {
        return false;
     }
     return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                                             & mpDefaultParameters->Numerical().KSPTolerances(),
                                             "KSPTolerances")->get().KSPRelative().present();
}

double HeartConfig::GetRelativeTolerance() const
{
	if (!GetUseRelativeTolerance())
    {
        EXCEPTION("Relative tolerance is not set in Chaste parameters");
    }   
    
	return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                           & mpDefaultParameters->Numerical().KSPTolerances(),
                           "KSPTolerances")->get().KSPRelative().get();
}

const char* HeartConfig::GetKSPSolver() const
{
    switch ( DecideLocation( & mpUserParameters->Numerical().KSPSolver(),
                             & mpDefaultParameters->Numerical().KSPSolver(),
                            "KSPSolver")->get() )
    {
    	case ksp_solver_type::gmres :
    		return "gmres";
		case ksp_solver_type::cg :
			return "cg";
		case ksp_solver_type::symmlq :
			return "symmlq";
    }
#define COVERAGE_IGNORE
    EXCEPTION("Unknown ksp solver");
#undef COVERAGE_IGNORE
}

const char* HeartConfig::GetKSPPreconditioner() const
{
    switch ( DecideLocation( & mpUserParameters->Numerical().KSPPreconditioner(),
                             & mpDefaultParameters->Numerical().KSPPreconditioner(),
                             "KSPPreconditioner")->get() )
    {
    	case ksp_preconditioner_type::ilu :
    		return "ilu";
    	case ksp_preconditioner_type::jacobi :
    		return "jacobi";
    	case ksp_preconditioner_type::bjacobi :
    		return "bjacobi";
    	
    }
#define COVERAGE_IGNORE
    EXCEPTION("Unknown ksp preconditioner");
#undef COVERAGE_IGNORE
}


/*
 *  Set methods
 */
// Simulation
void HeartConfig::SetSimulationDuration(double simulationDuration)
{
    mpUserParameters->Simulation().SimulationDuration().set(simulationDuration);
}

void HeartConfig::SetDomain(domain_type domain)
{
    mpUserParameters->Simulation().Domain().set(domain);
}

void HeartConfig::SetIonicModel(ionic_model_type ionicModel)
{
    mpUserParameters->Simulation().IonicModel().set(ionicModel);
}

void HeartConfig::SetMeshFileName(std::string meshPrefix, media_type fibreDefinition)
{
	if ( ! mpUserParameters->Simulation().Mesh().present())
	{
		mesh_type mesh_to_load;	
		mpUserParameters->Simulation().Mesh().set(mesh_to_load);	
	}
	
	mesh_type::LoadMesh::type mesh_prefix(meshPrefix, fibreDefinition);	
	mpUserParameters->Simulation().Mesh().get().LoadMesh().set(mesh_prefix);
}

void HeartConfig::SetOutputDirectory(std::string outputDirectory)
{
    mpUserParameters->Simulation().OutputDirectory().set(outputDirectory);
}

void HeartConfig::SetOutputFilenamePrefix(std::string outputFilenamePrefix)
{
    mpUserParameters->Simulation().OutputFilenamePrefix().set(outputFilenamePrefix);
}


// Physiological
void HeartConfig::SetIntracellularConductivities(const c_vector<double, 3>& intraConductivities)
{
    conductivities_type intra(intraConductivities[0],
                              intraConductivities[1],
                              intraConductivities[2]);

    mpUserParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 2>& intraConductivities)
{
    conductivities_type intra(intraConductivities[0],
                              intraConductivities[1],
                              DBL_MAX);

    mpUserParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 1>& intraConductivities)
{
    conductivities_type intra(intraConductivities[0],
                              DBL_MAX,
                              DBL_MAX);

    mpUserParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 3>& extraConductivities)
{
    conductivities_type extra(extraConductivities[0],
                              extraConductivities[1],
                              extraConductivities[2]);

    mpUserParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 2>& extraConductivities)
{
    conductivities_type extra(extraConductivities[0],
                              extraConductivities[1],
                              DBL_MAX);

    mpUserParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 1>& extraConductivities)
{
    conductivities_type extra(extraConductivities[0],
                              DBL_MAX,
                              DBL_MAX);

    mpUserParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetSurfaceAreaToVolumeRatio(double ratio)
{
    mpUserParameters->Physiological().SurfaceAreaToVolumeRatio().set(ratio);
}

void HeartConfig::SetCapacitance(double capacitance)
{
    mpUserParameters->Physiological().Capacitance().set(capacitance);
}


// Numerical
void HeartConfig::SetOdePdeAndPrintingTimeSteps(double odeTimeStep, double pdeTimeStep, double printingTimeStep)
{
    time_steps_type TimeSteps(odeTimeStep, pdeTimeStep, printingTimeStep);
    mpUserParameters->Numerical().TimeSteps().set(TimeSteps);
    CheckTimeSteps();
}

void HeartConfig::SetOdeTimeStep(double odeTimeStep)
{
    SetOdePdeAndPrintingTimeSteps(odeTimeStep, GetPdeTimeStep(), GetPrintingTimeStep());
}

void HeartConfig::SetPdeTimeStep(double pdeTimeStep)
{
    SetOdePdeAndPrintingTimeSteps(GetOdeTimeStep(), pdeTimeStep, GetPrintingTimeStep());
}

void HeartConfig::SetPrintingTimeStep(double printingTimeStep)
{
    SetOdePdeAndPrintingTimeSteps(GetOdeTimeStep(), GetPdeTimeStep(), printingTimeStep);
}

void HeartConfig::CheckTimeSteps() const
{
    if (GetOdeTimeStep() <= 0)
    {
        EXCEPTION("Ode time-step should be positive");
    }
    if (GetPdeTimeStep() <= 0)
    {
        EXCEPTION("Pde time-step should be positive");
    }
    if (GetPrintingTimeStep() <= 0.0)
    {
        EXCEPTION("Printing time-step should be positive");
    }
    if (GetPdeTimeStep()>GetPrintingTimeStep())
    {
        EXCEPTION("Printing time-step should not be smaller than PDE time step");
    }

    //If pde divides printing then the floating remainder ought to be close to
    //zero(+a smidge) or pde-a smidge
    double remainder=fmod(GetPrintingTimeStep(), GetPdeTimeStep());

    if ( remainder > DBL_EPSILON && remainder < GetPdeTimeStep()-DBL_EPSILON)
    {
        EXCEPTION("Printing time-step should a multiple of PDE time step");
    }
    
    if ( GetOdeTimeStep() > GetPdeTimeStep() )
    {
    	EXCEPTION("Ode time-step should not be greater than pde time-step");
    }
}


void HeartConfig::SetUseRelativeTolerance(double relativeTolerance)
{
    //Remove any reference to tolerances is user parameters
    mpUserParameters->Numerical().KSPTolerances().get().KSPAbsolute().reset();
    mpUserParameters->Numerical().KSPTolerances().get().KSPRelative().set(relativeTolerance);
}

void HeartConfig::SetUseAbsoluteTolerance(double absoluteTolerance)
{
    //Remove any reference to tolerances is user parameters
    mpUserParameters->Numerical().KSPTolerances().get().KSPRelative().reset();
    mpUserParameters->Numerical().KSPTolerances().get().KSPAbsolute().set(absoluteTolerance);
}

void HeartConfig::SetKSPSolver(const char* kspSolver)
{
	if ( strcmp(kspSolver, "gmres") == 0)
	{
		mpUserParameters->Numerical().KSPSolver().set(ksp_solver_type::gmres);
		return;
	}
	if ( strcmp(kspSolver, "cg") == 0)
	{
		mpUserParameters->Numerical().KSPSolver().set(ksp_solver_type::cg);
		return;
	}
	if ( strcmp(kspSolver, "symmlq") == 0)
	{
		mpUserParameters->Numerical().KSPSolver().set(ksp_solver_type::symmlq);
		return;
	}
	
	EXCEPTION("Unknown solver type provided");
}

void HeartConfig::SetKSPPreconditioner(const char* kspPreconditioner)
{
	if ( strcmp(kspPreconditioner, "ilu") == 0)
	{
        mpUserParameters->Numerical().KSPPreconditioner().set(ksp_preconditioner_type::ilu);
		return;
	}
	if ( strcmp(kspPreconditioner, "jacobi") == 0)
	{
        mpUserParameters->Numerical().KSPPreconditioner().set(ksp_preconditioner_type::jacobi);
		return;
	}
	if ( strcmp(kspPreconditioner, "bjacobi") == 0)
	{
        mpUserParameters->Numerical().KSPPreconditioner().set(ksp_preconditioner_type::bjacobi);
		return;
	}
	
	EXCEPTION("Unknown preconditioner type provided");
}
