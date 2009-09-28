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


#ifndef HEARTCONFIG_HPP_
#define HEARTCONFIG_HPP_

#include <string>
#include <vector>

#include "UblasIncludes.hpp"

#include "ArchiveLocationInfo.hpp"
#include "ChasteParameters.hpp"
#include "SimpleStimulus.hpp"
#include "ChasteCuboid.hpp"

#include <boost/shared_ptr.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>

// Needs to be included last
#include <boost/serialization/export.hpp>


/**
 * A singleton class containing configuration parameters for heart simulations.
 *
 * This class wraps the settings from the XML configuration file in a more friendly
 * interface, providing methods to read and write all the settings, and round-trip
 * them to/from XML format.  It also deals with the complexities of supporting
 * multiple versions of CodeSynthesis XSD.
 *
 * chaste_parameters_type is a convenience class created by CodeSynthesis XSD
 */
class HeartConfig
{
private:
    /**
     * Throws if the time steps don't obey constraints (within machine precision)
     * ode_step > 0.0
     * pde_step = n1 * ode_step (where n1 is a positive integer)
     * printing_step = n2 * pde_step (where n2 is a positive integer)
     */
    void CheckTimeSteps() const;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        //Only the Master should be writing the coonfiguration file
        if (PetscTools::AmMaster())
        {
            mpInstance->Write( true );
        }
        PetscTools::Barrier();
    }

    /**
     * Un-archive the object.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        std::string defaults_filename_xml = ArchiveLocationInfo::GetArchiveDirectory() + "ChasteDefaults.xml";
        HeartConfig::Instance()->SetDefaultsFile(defaults_filename_xml);

        std::string parameters_filename_xml = ArchiveLocationInfo::GetArchiveDirectory() + "ChasteParameters.xml";
        HeartConfig::Instance()->SetParametersFile(parameters_filename_xml);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()



public:
    /**
     * Call this method to access the global parameters holder.
     *
     * @return a single instance of the class
     */
    static HeartConfig* Instance();

    /**
     * @param useFixedSchemaLocation  whether to read the schema location from the XML
     *    file (false) or use the schema located at heart/src/io/ChasteParameters.xsd
     *    in the Chaste source tree (true).
     */
    void SetUseFixedSchemaLocation(bool useFixedSchemaLocation);

    /**
     * @param rFileName The name of the default file - set by default to "ChasteDefaults.xml" on construction
     */
    void SetDefaultsFile(const std::string& rFileName);
    /**
     * #mpUserParameters  is set to a new context associated with a parameters file
     * @param rFileName The name of the parameters file
     */
    void SetParametersFile(const std::string& rFileName);
    /**
     * Write out the complete configuration set (ChasteParameters
     * and ChasteDefaults) as an XML file.
     * Note that the location of ChasteParameters.xsd (schema definition)
     * will be hard-coded in the XML file.
     * @param useArchiveLocationInfo  if false, then use self's GetOutputDirectory() and open in /output subfolder
     *                                if true, then use ArchiveLocationInfo
     */
    void Write(bool useArchiveLocationInfo=false);

    /**
     * Throw away the current instance by resetting auto_ptr #mpInstance to NULL.
     * "New" another #mpInstance
     */
    static void Reset();

    /*
     *  Get methods
     */
     
    // Methods for asking the configuration file about which sections are defined. 
    /**
     *  Returns whether the configuration file defines a new simulation.
     * 
     *  @return is a new simulation? 
     */
    bool IsSimulationDefined() const;

    /**
     *  Returns whether the configuration file resumes an archived simulation.
     * 
     *  @return is a resumed simulation? 
     */
    bool IsSimulationResumed() const;         
     
    // Simulation
    unsigned GetSpaceDimension() const; /**< @return space dimension 1, 2 or 3.*/
    double GetSimulationDuration() const; /**< @return duration of the simulation (ms)*/
    /**
     * domain_type is an xsd convenience class type
     *
     * @return domain type of simulation bi- mono-domain
     */
    domain_type GetDomain() const;
    /**
     * Default cardiac cell model to use at all mesh nodes
     * (unless otherwise specified by IonicModelRegions)
     * ionic_models_available_type is an xsd convenience class type
     *
     * @return  type of model
     */
    ionic_models_available_type GetDefaultIonicModel() const;

    /**
     * Regions where we need to use a different cell model (think infarction)
     * ionic_models_available_type is an xsd convenience class type.
     *
     * \todo - do we assume the vectors are initially empty?
     * The standard vectors returned are of the same length (one entry per region)
     *
     * @param definedRegions vector of axis-aligned box regions (one per cellular heterogeneity)
     * @param ionicModels vector of models (one per cellular heterogeneity)
     * \todo No set method
     */
     void GetIonicModelRegions(std::vector<ChasteCuboid>& definedRegions,
                               std::vector<ionic_models_available_type>& ionicModels) const;


    bool IsMeshProvided() const; /**< @return true if a mesh file name is given.  (Otherwise it's assumed that this is a cuboid simulation.)*/
    bool GetCreateMesh() const; /**< @return true if it's cuboid simulation (no mesh on disk)*/
    bool GetCreateSlab() const; /**< @return true if it's cuboid simulation (no mesh on disk)*/
    bool GetCreateSheet() const; /**< @return true if it's cuboid simulation (no mesh on disk)*/
    bool GetCreateFibre() const; /**< @return true if it's cuboid simulation (no mesh on disk)*/
    bool GetLoadMesh() const; /**< @return true if a mesh file name is given and we are expecting to load a mesh from file*/
    ///\todo IsMeshProvided and GetLoadMesh are subtly different but very similar.  Can one of them go?
    /**
     * @param slabDimensions  return vector for the (cuboid) mesh dimensions (cm)
     */
    void GetSlabDimensions(c_vector<double, 3>& slabDimensions) const;
    /**
     * @param sheetDimensions  return vector for the (cuboid) mesh dimensions (cm)
     */
    void GetSheetDimensions(c_vector<double, 2>& sheetDimensions) const;
    /**
     * @param fibreLength  return vector for the (cuboid) mesh dimensions (cm)
     */
    void GetFibreLength(c_vector<double, 1>& fibreLength) const;
    double GetInterNodeSpace() const; /**< @return internode space of cuboid mesh (cm)*/

    std::string GetMeshName() const;/**< @return path/basename of mesh files*/

    media_type GetConductivityMedia() const;/**< @return media (Orthotropic/Axisymmetric/NoFibreOrientation) so that we know whether to read a .ortho/.axi file*/

    /**
     * Return a number of stimulated regions (Axis-aligned boxes)
     * \todo - do we assume the vectors are initially empty?
     * The returned std::vectors are all of the same length
     * @param rStimuliApplied  rStimuliApplied[0] is stimulus for the first region
     * @param rStimulatedAreas  rStimulatedAreas[0] is the first region to be stimulated
     *
     * \todo There is no set method
     */
    void GetStimuli(std::vector<boost::shared_ptr<SimpleStimulus> >& rStimuliApplied, std::vector<ChasteCuboid>& rStimulatedAreas) const;

    /**
     * Return a number of heterogeneous regions (Axis-aligned boxes) for special gating variable changes
     * \todo - do we assume the vectors are initially empty?
     * The returned std::vectors are all of the same length
     * @param cellHeterogeneityAreas  cellHeterogeneityAreas[0] is the first region
     * @param scaleFactorGks  scaleFactorGks[0] is a scaling factorfor the first region
     * @param scaleFactorIto  scaleFactorIto[0] is a scaling factorfor the first region
     * @param scaleFactorGkr  scaleFactorGkr[0] is a scaling factorfor the first region
     * \todo There is no set method
     */
    void GetCellHeterogeneities(std::vector<ChasteCuboid>& cellHeterogeneityAreas,
                                std::vector<double>& scaleFactorGks,
                                std::vector<double>& scaleFactorIto,
                                std::vector<double>& scaleFactorGkr) const;
    bool GetConductivityHeterogeneitiesProvided() const; /**< @return  true if there are conductivity heterogeneities for GetConductivityHeterogeneities to return*/
    /**
     * Return a number of heterogeneous regions (Axis-aligned boxes)
     * \todo - do we assume the vectors are initially empty?
     * The returned std::vectors are all of the same length
     * @param conductivitiesHeterogeneityAreas  conductivitiesHeterogeneityAreas[0] is the first region
     * @param intraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param extraConductivities  extraConductivities[0] is conductivity vector for the first region
     */
    void GetConductivityHeterogeneities(std::vector<ChasteCuboid>& conductivitiesHeterogeneityAreas,
                                        std::vector< c_vector<double,3> >& intraConductivities,
                                        std::vector< c_vector<double,3> >& extraConductivities) const;
    std::string GetOutputDirectory() const; /**< @return output directory path name*/

    /**
     * @return  Prefix for files
     * If set to "res" this produces
     * [path]/res.h5
     * [path]/output/res_mesh.pts
     * [path]/output/res_mesh.tri
     * [path]/output/res_parameters.xml  (a copy of this configuration at the end of the simulation)
     * [path]/output/res_times.info
     * [path]/output/res_V.dat
     */
    std::string GetOutputFilenamePrefix() const;

    /**
     * @return true any extra output variables have been requested
     */
    bool GetOutputVariablesProvided() const;

    /**
     * Get the extra output variables from the xml file
     *
     * @param outputVariables reference to std::vector to contain the output variables requested
     */
    void GetOutputVariables(std::vector<std::string> &outputVariables) const;

    /**
     * Get whether simulation should be archived or not
     * 
     * @return archive simulation
     */
    bool GetSaveSimulation() const;


    // Physiological
    /**
     * 3D version
     * @param intraConductivities  DIM-vector for returning intracellular conductivities (mS/cm)
     */
    void GetIntracellularConductivities(c_vector<double, 3>& intraConductivities) const;
    /**
     * 2D version
     * @param intraConductivities  DIM-vector for returning intracellular conductivities (mS/cm)
     */
    void GetIntracellularConductivities(c_vector<double, 2>& intraConductivities) const;
    /**
     * 1D version
     * @param intraConductivities  DIM-vector for returning intracellular conductivities (mS/cm)
     */
    void GetIntracellularConductivities(c_vector<double, 1>& intraConductivities) const;

    /**
     * 3D version
     * @param extraConductivities  DIM-vector for returning extracellular conductivities (mS/cm)
     */
    void GetExtracellularConductivities(c_vector<double, 3>& extraConductivities) const;
    /**
     * 2D version
     * @param extraConductivities  DIM-vector for returning extracellular conductivities (mS/cm)
     */
    void GetExtracellularConductivities(c_vector<double, 2>& extraConductivities) const;
    /**
     * 1D version
     * @param extraConductivities  DIM-vector for returning extracellular conductivities (mS/cm)
     */
    void GetExtracellularConductivities(c_vector<double, 1>& extraConductivities) const;

    double GetBathConductivity() const; /**< @return conductivity for perfusing bath (mS/cm)*/


    double GetSurfaceAreaToVolumeRatio() const; /**< @return surface area to volume ratio chi a.k.a Am for PDE (1/cm)*/

    double GetCapacitance() const; /**< @return surface capacitance Cm for PDE (uF/cm^2)*/

    // Numerical
    double GetOdeTimeStep() const; /**< @return ODE time-step (ms)*/
    double GetPdeTimeStep() const; /**< @return PDE time-step (ms)*/
    double GetPrintingTimeStep() const; /**< @return priting time-step (ms)*/

    bool GetUseAbsoluteTolerance() const; /**< @return true if we are using KSP absolute tolerance*/
    double GetAbsoluteTolerance() const; /**< @return KSP absolute tolerance (or throw if we are using relative)*/

    bool GetUseRelativeTolerance() const; /**< @return true if we are using KSP relative tolerance*/
    double GetRelativeTolerance() const;  /**< @return KSP relative tolerance (or throw if we are using absolute)*/

    const char* GetKSPSolver() const; /**< @return name of -ksp_type from {"gmres", "cg", "symmlq"}*/
    const char* GetKSPPreconditioner() const; /**< @return name of -pc_type from {"ilu", "jacobi", "bjacobi", "hypre", "none"}*/

    // Post processing
    /**
     * @return true if there is a post-processing section
     */
    bool IsPostProcessingSectionPresent() const;

    /**
     * @return true if any post-processing information has been requested
     */
    bool IsPostProcessingRequested() const;
    
    /**
     * @return true if APD maps have been requested
     */
    bool IsApdMapsRequested() const;
    /**
     * @param apdMaps  each entry is a request for a map with
     *  - a threshold (in mV)
     *  - a percentage in the range [1, 100)
     */
    void GetApdMaps(std::vector<std::pair<double,double> >& apdMaps) const;

    /**
     * @return true if upstroke time maps have been requested
     */
    bool IsUpstrokeTimeMapsRequested() const;
    /**
     * @param upstrokeTimeMaps  each entry is a request for a map with
     *  - a threshold (in mV)
     */
    void GetUpstrokeTimeMaps (std::vector<double>& upstrokeTimeMaps) const;

    /**
     * @return true maximum upstroke velocity maps have been requested
     */
    bool IsMaxUpstrokeVelocityMapRequested() const;

    /**
     * @param upstrokeVelocityMaps  each entry is a request for a map with
     *  - a threshold (in mV, defaulted to -30 mV)
     */
    void GetMaxUpstrokeVelocityMaps(std::vector<double>& upstrokeVelocityMaps) const;

    /**
     * @return true if conduction velocity maps have been requested
     */
    bool IsConductionVelocityMapsRequested() const;

    /**
     * @param conductionVelocityMaps  each entry is a request for a map with
     *  - an index to treat as ths source for wave propagation
     */
    void GetConductionVelocityMaps(std::vector<unsigned>& conductionVelocityMaps) const;


     /*
     *  Set methods
     */

    // Simulation
    /** Set the configuration dimension
     * @param spaceDimension 1, 2 or 3.
     */
    void SetSpaceDimension(unsigned spaceDimension);

    /** Set the configuration simulation duration
     * @param simulationDuration duration of the simulation (ms)
     */
    void SetSimulationDuration(double simulationDuration);
    /**
     * Set the configuration to run mono or bidomain
     * domain_type is an xsd convenience class type
     *
     * @param domain type of simulation bi- mono-domain
     */
    void SetDomain(domain_type domain);
    /**
     * Set the configuration to place the given cardiac cell models at all mesh nodes
     * (unless otherwise specified by IonicModelRegions)
     * ionic_models_available_type is an xsd convenience class type
     *
     * @param ionicModel  type of model
     */
    void SetDefaultIonicModel(ionic_models_available_type ionicModel);

    /**
     * Set dimensions of simulation for use with a cuboid mesh generated on the fly.  3-D.
     * @param x  length in 1st dimension (cm)
     * @param y  length in 2nd dimension (cm)
     * @param z  length in 3rd dimension (cm)
     * @param inter_node_space  Spacing in cartesian direction (cm). Diagonals will be longer.
     */
    void SetSlabDimensions(double x, double y, double z, double inter_node_space);
    /**
     * Set dimensions of simulation for use with a cuboid mesh generated on the fly.  2-D.
     * @param x  length in 1st dimension (cm)
     * @param y  length in 2nd dimension (cm)
     * @param inter_node_space  Spacing in cartesian direction (cm). Diagonals will be longer.
     */
    void SetSheetDimensions(double x, double y, double inter_node_space);
    /**
     * Set dimensions of simulation for use with a cuboid mesh generated on the fly.  1-D.
     * @param x  length in 1st dimension (cm)
     * @param inter_node_space  Spacing in cartesian direction (cm).
     */
    void SetFibreLength(double x, double inter_node_space);

    /**
     * Sets the name of a mesh to be read from disk for this simulation
     * @param meshPrefix  path and basename of a set of mesh files (.nodes .ele etc) in triangle/tetget format
     * @param fibreDefinition  if set (Orthotropic/Axisymmetric) then a (.ortho/.axi) file should also be read
     * \todo There is no Get method
     */
    void SetMeshFileName(std::string meshPrefix, media_type fibreDefinition=media_type::NoFibreOrientation);

    /**
     * Set a number of heterogeneous regions (Axis-aligned boxes)
     * It is assumed that the std::vectors are all of the same length
     * @param cornerA  cornerA[0] is the lowest vertex of the first region
     * @param cornerB  cornerB[0] is the highest vertex of the first region
     * @param intraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param extraConductivities  extraConductivities[0] is conductivity vector for the first region
     */
    void SetConductivityHeterogeneities(std::vector< c_vector<double,3> >& cornerA,
                                        std::vector< c_vector<double,3> >& cornerB,
                                        std::vector< c_vector<double,3> >& intraConductivities,
                                        std::vector< c_vector<double,3> >& extraConductivities);

    /**
     * @param rOutputDirectory  Full path to output directory (will be created if necessary)
     */
    void SetOutputDirectory(const std::string& rOutputDirectory);
    /**
     * @param rOutputFilenamePrefix  Prefix for files
     * If set to "res" this will produce
     * [path]/res.h5
     * [path]/output/res_mesh.pts
     * [path]/output/res_mesh.tri
     * [path]/output/res_parameters.xml  (a copy of this configuration at the end of the simulation)
     * [path]/output/res_times.info
     * [path]/output/res_V.dat
     */
    void SetOutputFilenamePrefix(const std::string& rOutputFilenamePrefix);

    /**
     * @param rOutputVariables  a vector of std::strings of the names
     * of each variable that should be outputted at each time step.
     *
     * USING THIS METHOD WILL OVERRIDE THE ANY OUTPUT VARIABLES SET IN THE XML FILE
     */
    void SetOutputVariables(const std::vector<std::string>& rOutputVariables);

    /**
     * Set whether the simulation should be archived or not
     * 
     * @param saveSimulation archive simulation
     */
     void SetSaveSimulation(bool saveSimulation);


    // Physiological
    /**
     * 3D version
     * @param intraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 3>& intraConductivities);
    /**
     * 2D version
     * @param intraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 2>& intraConductivities);
    /**
     * 1D version
     * @param intraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 1>& intraConductivities);

    /**
     * 3D version
     * @param extraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 3>& extraConductivities);
    /**
     * 2D version
     * @param extraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 2>& extraConductivities);
    /**
     * 1D version
     * @param extraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 1>& extraConductivities);

    /**
     * Set bath conductivity
     * @param bathConductivity conductivity for perfusing bath (mS/cm)
     * \todo Is this used anywhere?
     */
    void SetBathConductivity(double bathConductivity);

    /**
     * Set surface area to volume ratio Am (for PDE)
     * @param ratio (1/cm)
     */
    void SetSurfaceAreaToVolumeRatio(double ratio);

    /**
     * Set surface capacitance Cm (for PDE)
     * @param capacitance (uF/cm^2)
     */
    void SetCapacitance(double capacitance);

    // Numerical
    /** Set the configuration to use ode, pde and printing times of given values
     * Calls CheckTimeSteps to ensure compatibility
     * @param odeTimeStep  ode value to use
     * @param pdeTimeStep  pde value to use
     * @param printingTimeStep  printing value to use
     */
    void SetOdePdeAndPrintingTimeSteps(double odeTimeStep, double pdeTimeStep, double printingTimeStep);
   /** Set the configuration to use ode time of given value
     * Calls CheckTimeSteps via SetOdePdeAndPrintingTimeSteps
     * @param odeTimeStep  the value to use
     */
    void SetOdeTimeStep(double odeTimeStep);
    /** Set the configuration to use pde time of given value
     * Calls CheckTimeSteps via SetOdePdeAndPrintingTimeSteps
     * @param pdeTimeStep  the value to use
     */
    void SetPdeTimeStep(double pdeTimeStep);

    /** Set the configuration to use printing time of given value
     * Calls CheckTimeSteps via SetOdePdeAndPrintingTimeSteps
     * @param printingTimeStep  the value to use
     */
     void SetPrintingTimeStep(double printingTimeStep);

    /** Set the configuration to use KSP relative tolerance of given value
     * @param relativeTolerance  the value to use
     */
    void SetUseRelativeTolerance(double relativeTolerance);
    /** Set the configuration to use KSP absolute tolerance of given value
     * @param absoluteTolerance  the value to use
     */
    void SetUseAbsoluteTolerance(double absoluteTolerance);

    /** Set the type of KSP solver as with the flag "-ksp_type"
     * @param kspSolver  a string from {"gmres", "cg", "symmlq"}
     */
    void SetKSPSolver(const char* kspSolver);
    /** Set the type of preconditioner as with the flag "-pc_type"
     * @param kspPreconditioner  a string from {"ilu", "jacobi", "bjacobi", "hypre", "blockdiagonal", "none"}
     */
    void SetKSPPreconditioner(const char* kspPreconditioner);

    /** Set the parameters of the apd map requested
     *
     *  @param apdMaps  each entry is a request for a map with
     *  - a threshold (in mV)
     *  - a percentage in the range [1, 100) (ranges are not checked by this method, but during the calculation)
     */
    void SetApdMaps(const std::vector<std::pair<double,double> >& apdMaps);

    /** Set the parameters of the upstroke time map requested
     *
     *  @param upstrokeTimeMaps  is the list of thresholds (? ///\todo improve the description of threshold) with respect to which the upstroke time maps are calculated.
     */
    void SetUpstrokeTimeMaps (std::vector<double>& upstrokeTimeMaps);

    /** Set the parameters of the maximal upstroke velocity map requested
     *
     *  @param maxUpstrokeVelocityMaps is the list of thresholds (? ///\todo improve the description of threshold) with respect to which the upstroke velocity maps are calculated. 
     */
    void SetMaxUpstrokeVelocityMaps (std::vector<double>& maxUpstrokeVelocityMaps);

    /** Set the parameters of the conduction velocity map requested
     *
     *  @param conductionVelocityMaps is a list of origin node indices. One map is created for each origin node.
     */
    void SetConductionVelocityMaps (std::vector<unsigned>& conductionVelocityMaps);

    ~HeartConfig(); /**< Destructor*/
protected:
    // Only to be accessed by the tests
    friend class TestHeartConfig;

private:
    /*Constructor is private, since the class is only accessed by the singleton instance() method*/
    HeartConfig();

    /** Pointer to parameters read from the user's input XML file
     * (override those given by #mpDefaultParameters).
     */
    boost::shared_ptr<chaste_parameters_type> mpUserParameters;
    /** Pointer to parameters read from the default input XML file (to be read before
     * #mpUserParameters, but may be subsequently overridden).
     */
    boost::shared_ptr<chaste_parameters_type> mpDefaultParameters;

    /** The single instance of the class */
    static std::auto_ptr<HeartConfig> mpInstance;

    /**
     * Whether to read the schema location from the XML file (false) or use the schema
     * located at heart/src/io/ChasteParameters.xsd in the Chaste source tree (true).
     */
    bool mUseFixedSchemaLocation;

    /**
     * DecideLocation is a convenience method used to get the correct parameter value
     * from the defaults/parameters files.  It checks if the first value  is present and (if not)
     * moves onto the second
     *
     * @param params_ptr  Pointer to quantity within the parameters file (checked first, since it will override a default)
     * @param defaults_ptr  Pointer to quantity within the defaults file (used if there was no override)
     * @param nameParameter Name of quatity within params_ptr/defaults_ptr (so we can throw a meaningful exception if it's not found)
     */
    template<class TYPE>
    TYPE* DecideLocation(TYPE* params_ptr, TYPE* defaults_ptr, const std::string& nameParameter) const;

    /**
     * CheckSimulationIsDefined is a convience method for checking if the "<"Simulation">" element
     * has been defined and therefore is safe to use the Simulation().get() pointer to access
     * other data.
     * 
     * Throws and exception if not.
     * 
     * @param callingMethod string describing the get method performing the check.
     */
     void CheckSimulationIsDefined(std::string callingMethod="") const;

    /**
     * Utility method to parse an XML parameters file.
     * @param rFileName  Name of XML file
     */
    boost::shared_ptr<chaste_parameters_type> ReadFile(const std::string& rFileName);

};

#include "TemplatedExport.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(HeartConfig);

#endif /*HEARTCONFIG_HPP_*/
