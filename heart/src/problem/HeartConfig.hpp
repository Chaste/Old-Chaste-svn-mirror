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


#ifndef HEARTCONFIG_HPP_
#define HEARTCONFIG_HPP_

#include <string>
#include <vector>

#include "UblasIncludes.hpp"

#include "ArchiveLocationInfo.hpp"
#include "ChasteParameters_2_1.hpp"

#include "AbstractStimulusFunction.hpp"
// These are needed here for Boost < 1.37
#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"

#include "ChasteCuboid.hpp"
#include "ChasteEllipsoid.hpp"
#include "AbstractTetrahedralMesh.hpp"

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>


#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>

namespace cp = chaste::parameters::v2_1;

// Forward declaration to avoid circular includes
class HeartFileFinder;


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
        //Only the Master should be writing the configuration file
        if (PetscTools::AmMaster())
        {
            mpInstance->Write( true );
        }
        PetscTools::Barrier("HeartConfig::save");
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
        /*
         *  This method implements the logic required by HeartConfig to be able to handle resuming a simulation via the executable.
         *
         *  When the control reaches the method mpUserParameters and mpDefaultParameters point to the files specified as resuming parameters.
         *  However SetDefaultsFile() and SetParametersFile() will set those variables to point to the archived parameters.
         *
         *  We make a temporary copy of mpUserParameters so we don't lose its content. At the end of the method we update the new mpUserParameters
         *  with the resuming parameters.
         */
        assert(mpUserParameters.use_count() > 0);
        boost::shared_ptr<cp::chaste_parameters_type> p_new_parameters = mpUserParameters;

        std::string defaults_filename_xml = ArchiveLocationInfo::GetArchiveDirectory() + "ChasteDefaults.xml";
        HeartConfig::Instance()->SetDefaultsFile(defaults_filename_xml);

        /*
         *  When we unarchive a simulation, we load the old parameters file in order to inherit things such
         *  as default cell model, stimuli, heterogeneities, ... This has the side effect of inheriting the
         *  <CheckpointSimulation> element (if defined).
         *
         *  We disable checkpointing definition coming from the unarchived config file. We will enable it again
         *  if defined in the resume config file.
         */
        std::string parameters_filename_xml = ArchiveLocationInfo::GetArchiveDirectory() + "ChasteParameters.xml";
        HeartConfig::Instance()->SetParametersFile(parameters_filename_xml);

        HeartConfig::Instance()->SetCheckpointSimulation(false);

        // If we are resuming a simulation, some parameters can be altered at this point.
        if (p_new_parameters->ResumeSimulation().present())
        {
            UpdateParametersFromResumeSimulation(p_new_parameters);
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
    /**
     * When loading a simulation from archive, some parameters can get overridden by the content of the ResumeSimulation
     * element.  This method does that.
     * 
     * @param pResumeParameters  the parameters containing the ResumeSimulation element.
     */
    void UpdateParametersFromResumeSimulation(boost::shared_ptr<cp::chaste_parameters_type> pResumeParameters);

public:
    /**
     * Our type for specifying schema location properties: a map from namespace URI
     * to schema URI.  The default namespace is specified by an empty namespace URI.
     */
    typedef std::map<std::string, std::string> SchemaLocationsMap;

private:
    /**
     * Fixed location of schema files for the different Chaste parameters namespaces.
     */
    SchemaLocationsMap mSchemaLocations;

    /**
     * Set default schema locations in the Chaste source tree.
     */
    void SetDefaultSchemaLocations();

    /**
     * Helper method for URL-escaping spaces in file paths, to avoid confusing Xerces
     * regarding schema locations.  Note that this is a very specific fix: it doesn't
     * do general URL-escaping.
     *
     * @param rPath  the path to escape
     */
    std::string EscapeSpaces(const std::string& rPath);

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
     *    in the Chaste source tree (or specified with SetFixedSchemaLocations()) (true).
     */
    void SetUseFixedSchemaLocation(bool useFixedSchemaLocation);

    /**
     * Set the schema files to use.
     * Also calls SetUseFixedSchemaLocation(true).
     *
     * @param rSchemaLocations  map from namespace URI to schema URI
     */
    void SetFixedSchemaLocations(const SchemaLocationsMap& rSchemaLocations);

    /**
     * Allows users to override the built-in defaults.
     * @param rFileName The name of the defaults file
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
     * @param useArchiveLocationInfo  if false, then use self's GetOutputDirectory() and open in *named* subfolder
     *                                if true, then use ArchiveLocationInfo
     * @param subfolderName -- where to store with respect to GetOutputDirectory()
     *
     * @note This method is collective if useArchiveLocationInfo is false
     */
    void Write(bool useArchiveLocationInfo=false, std::string subfolderName="output");

    /**
     * Try to copy the latest version of the schema to the given directory.
     * If we can't find the latest version of the schema, generate a warning.
     * @param rToDirectory  directory to copy to
     */
    void CopySchema(const std::string& rToDirectory);

    /**
     * Utility method to parse an XML parameters file.
     * @param rFileName  Name of XML file
     */
    boost::shared_ptr<cp::chaste_parameters_type> ReadFile(const std::string& rFileName);

    /**
     * Throw away the current instance by resetting auto_ptr #mpInstance to NULL.
     * "New" another #mpInstance
     */
    static void Reset();
    
    ~HeartConfig(); /**< Destructor*/
    
    /**
     * Get the Chaste version of a parameters file, given its namespace URI.
     * The version will be encoded as major*1000+minor.
     * 
     * @param rNamespaceUri  the namespace URI of the parameters file
     */
    unsigned GetVersionFromNamespace(const std::string& rNamespaceUri);

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
     * cp::domain_type is an xsd convenience class type
     *
     * @return domain type of simulation: bi- or mono-domain
     */
    cp::domain_type GetDomain() const;

    /**
     * Default cardiac cell model to use at all mesh nodes
     * (unless otherwise specified by GetIonicModelRegions).
     * cp::ionic_model_selection_type is generated automatically from the XML Schema.
     *
     * @return  type of model
     */
    cp::ionic_model_selection_type GetDefaultIonicModel() const;

    /**
     * Regions where we need to use a different cell model (think infarction).
     * cp::ionic_model_selection_type is generated automatically from the XML Schema.
     *
     * The supplied vectors are first cleared, then filled in with the information from the
     * parameters files.  On return, both vectors will be the same length (one entry per region).
     *
     * @param rDefinedRegions vector of axis-aligned box regions (one per cellular heterogeneity)
     * @param rIonicModels vector of models (one per cellular heterogeneity)
     */
     template<unsigned DIM>
     void GetIonicModelRegions(std::vector<ChasteCuboid<DIM> >& rDefinedRegions,
                               std::vector<cp::ionic_model_selection_type>& rIonicModels) const;

    /**
     * Set the regions where we need to use a different cell model (think infarction).
     * Unlike the get method, this is currently only supported in 3d.
     * cp::ionic_model_selection_type is generated automatically from the XML Schema.
     *
     * The input standard vectors must be of the same length (one entry per region)
     * otherwise the method throws.
     *
     * @param rDefinedRegions vector of axis-aligned box regions (one per cellular heterogeneity)
     * @param rIonicModels vector of models (one per cellular heterogeneity)
     */
     void SetIonicModelRegions(std::vector<ChasteCuboid<3> >& rDefinedRegions,
                               std::vector<cp::ionic_model_selection_type>& rIonicModels) const;

    bool IsMeshProvided() const; /**< @return true if a mesh file name is given.  (Otherwise it's assumed that this is a cuboid simulation.)*/
    bool GetCreateMesh() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
    bool GetCreateSlab() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
    bool GetCreateSheet() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
    bool GetCreateFibre() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
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

    cp::media_type GetConductivityMedia() const;/**< @return media (Orthotropic/Axisymmetric/NoFibreOrientation) so that we know whether to read a .ortho/.axi file*/

    /**
     * Return a number of stimulated regions (Axis-aligned boxes)
     * \todo - do we assume the vectors are initially empty?
     * The returned std::vectors are all of the same length
     * @param rStimuliApplied  rStimuliApplied[0] is stimulus for the first region
     * @param rStimulatedAreas  rStimulatedAreas[0] is the first region to be stimulated
     *
     * \todo There is no set method
     */
     template<unsigned DIM>
    void GetStimuli(std::vector<boost::shared_ptr<AbstractStimulusFunction> >& rStimuliApplied, std::vector<ChasteCuboid<DIM> >& rStimulatedAreas) const;

    /**
     * Reads from the XML file the cellular hetrogeneities. It fugures out whether the user specified a cuboid
     * or a transmural-type of hetrogeneities. In the latter case, it stores the percentage values of Epi and Endo layers
     * in two member variables, accessible via get methods. It also checks if the user-supplied numbers are consistent (i.e., positive and add up to less than 1)
     * Return a number of heterogeneous regions for special gating variable changes
     * \todo - do we assume the vectors are initially empty?
     * The returned std::vectors are all of the same length
     * @param rCellHeterogeneityRegions  cellHeterogeneityAreas[0] is the first region
     * @param rScaleFactorGks  scaleFactorGks[0] is a scaling factor for the first region
     * @param rScaleFactorIto  scaleFactorIto[0] is a scaling factor for the first region
     * @param rScaleFactorGkr  scaleFactorGkr[0] is a scaling factor for the first region
     * @param pParameterSettings  specification of named parameters to set on the cell models; each entry is a map
     *     from parameter name to value.
     * \todo There is no set method
     */
    template<unsigned DIM>
    void GetCellHeterogeneities( std::vector<AbstractChasteRegion<DIM>* >& rCellHeterogeneityRegions,
                                 std::vector<double>& rScaleFactorGks,
                                 std::vector<double>& rScaleFactorIto,
                                 std::vector<double>& rScaleFactorGkr,
                                 std::vector<std::map<std::string, double> >* pParameterSettings);

    bool GetConductivityHeterogeneitiesProvided() const; /**< @return  true if there are conductivity heterogeneities for GetConductivityHeterogeneities to return*/

    /**
     * @return the value of the flag that tells whether the user asked for cellular transmural heterogeneities
     */
    bool AreCellularTransmuralHeterogeneitiesRequested();

    /**
     * @return the value of the flag that tells whether the user asked for cellular heterogeneities with cuboids
     */
    bool AreCellularHeterogeneitiesSpecifiedByCuboids();

    /**
     * @return the fraction of epicardial layer
     */
    double GetEpiLayerFraction();

    /**
     * @return the fraction of endocardial layer
     */
    double GetEndoLayerFraction();

    /**
     * @return the fraction of endocardial layer
     */
    double GetMidLayerFraction();

    /**
     * @return the index with which the epicardial layer is supplied (i.e., the order it comes in the XML file)
     */
    unsigned GetEpiLayerIndex();

    /**
     * @return the index with which the endocardial layer is supplied (i.e., the order it comes in the XML file)
     */
    unsigned GetEndoLayerIndex();

    /**
     * @return the index with which the midmyocardial layer is supplied (i.e., the order it comes in the XML file)
     */
    unsigned GetMidLayerIndex();


    /**
     * Return a number of heterogeneous regions (Axis-aligned boxes)
     * \todo - do we assume the vectors are initially empty?
     * The returned std::vectors are all of the same length
     * @param conductivitiesHeterogeneityAreas  conductivitiesHeterogeneityAreas[0] is the first region
     * @param intraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param extraConductivities  extraConductivities[0] is conductivity vector for the first region
     */
    template<unsigned DIM>
    void GetConductivityHeterogeneities(std::vector<AbstractChasteRegion<DIM>* >& conductivitiesHeterogeneityAreas,
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
     * @return true iff any extra output variables have been requested
     */
    bool GetOutputVariablesProvided() const;

    /**
     * Get the extra output variables from the xml file.
     *
     * @param rOutputVariables reference to std::vector to contain the output variables requested.
     *    Note: will be cleared before being filled.
     */
    void GetOutputVariables(std::vector<std::string>& rOutputVariables) const;
    
    /**
     * @return whether to write output HDF5 file using the original
     * mesh permutation (in situations where a parallel partition may have
     * permuted the node).  The default is to use the new, not original permutation,
     */
    bool GetOutputUsingOriginalNodeOrdering();

    /**
     * Get whether simulation should be checkpointed or not
     *
     * @return archive simulation
     */
    bool GetCheckpointSimulation() const;

    /**
     * Get checkpointing timestep
     *
     * @return checkpointing timestep
     */
    double GetCheckpointTimestep() const;

    /**
     * Get number of checkpoints to keep on disk
     *
     * @return checkpointing timestep
     */
    unsigned GetMaxCheckpointsOnDisk() const;

    // ResumeSimulation
    /**
     * Get directory where the archived simulation to resume is defined
     */
    HeartFileFinder GetArchivedSimulationDir() const;


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
    const char* GetKSPPreconditioner() const; /**< @return name of -pc_type from {"jacobi", "bjacobi", "hypre", "ml", "spai", "blockdiagonal", "ldufactorisation", "none"}*/


    // Adaptivity
    /**
     * @return true if there is an adaptivity section
     */
    bool IsAdaptivityParametersPresent() const;

    /**
     * @return value of target error used in adaptivity
     */
    double GetTargetErrorForAdaptivity() const;

    /**
     * @return value of sigma used in adaptivity
     */
    double GetSigmaForAdaptivity() const;

    /**
     * @return maximum edge length in mesh after an adapt
     */
    double GetMaxEdgeLengthForAdaptivity() const;

    /**
     * @return minimum edge length in mesh after an adapt
     */
    double GetMinEdgeLengthForAdaptivity() const;

    /**
     * @return value of gradation used in adaptivity
     */
    double GetGradationForAdaptivity() const;

    /**
     * @return maximum number of nodes in mesh after an adapt
     */
    unsigned GetMaxNodesForAdaptivity() const;

    /**
     * @return number of adaptive sweeps through mesh during adaptivity
     */
    unsigned GetNumberOfAdaptiveSweeps() const;

    // Post processing
    /**
     * @return true if there is a post-processing section
     */
    bool IsPostProcessingSectionPresent() const;
    
    /**
     * Create a PostProcessing section in the user parameters if one doesn't exist.
     */
    void EnsurePostProcessingSectionPresent();

    /**
     * @return true if any post-processing information has been requested
     */
    bool IsPostProcessingRequested() const;

    /**
     * @return true if APD maps have been requested
     */
    bool IsApdMapsRequested() const;
    /**
     * @param rApdMaps  each entry is a request for a map with
     *  - a percentage in the range [1, 100)
     *  - a threshold (in mV)
     */
    void GetApdMaps(std::vector<std::pair<double,double> >& rApdMaps) const;

    /**
     * @return true if upstroke time maps have been requested
     */
    bool IsUpstrokeTimeMapsRequested() const;
    /**
     * @param rUpstrokeTimeMaps  each entry is a request for a map with
     *  - a threshold (in mV)
     */
    void GetUpstrokeTimeMaps (std::vector<double>& rUpstrokeTimeMaps) const;

    /**
     * @return true maximum upstroke velocity maps have been requested
     */
    bool IsMaxUpstrokeVelocityMapRequested() const;

    /**
     * @param rUpstrokeVelocityMaps  each entry is a request for a map with
     *  - a threshold (in mV, defaulted to -30 mV)
     */
    void GetMaxUpstrokeVelocityMaps(std::vector<double>& rUpstrokeVelocityMaps) const;

    /**
     * @return true if conduction velocity maps have been requested
     */
    bool IsConductionVelocityMapsRequested() const;

    /**
     * @param rConductionVelocityMaps  each entry is a request for a map with
     *  - an index to treat as the source for wave propagation
     */
    void GetConductionVelocityMaps(std::vector<unsigned>& rConductionVelocityMaps) const;


    // Output visualization

    /** Whether there is an OutputVisualizer element present. */
    bool IsOutputVisualizerPresent() const;

    /** Whether to convert the output from HDF5 to meshalyzer readable format */
    bool GetVisualizeWithMeshalyzer() const;

    /** Whether to convert the output from HDF5 to Cmgui readable format */
    bool GetVisualizeWithCmgui() const;

    /** Whether to convert the output from HDF5 to Vtk readable format */
    bool GetVisualizeWithVtk() const;
    
    /**
     * @return true if there is an electrodes section
     */
    bool IsElectrodesPresent() const;    

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
     * cp::domain_type is an xsd convenience class type
     *
     * @param rDomain type of simulation bi- mono-domain
     */
    void SetDomain(const cp::domain_type& rDomain);

    /**
     * Set the configuration to place the given cardiac cell models at all mesh nodes
     * (unless otherwise specified by SetIonicModelRegions).
     * cp::ionic_models_available_type is generated automatically from the XML Schema.
     *
     * @param rIonicModel  type of model
     */
    void SetDefaultIonicModel(const cp::ionic_models_available_type& rIonicModel);

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
    void SetMeshFileName(std::string meshPrefix, cp::media_type fibreDefinition=cp::media_type::NoFibreOrientation);

    /**
     * Set a number of heterogeneous regions (Axis-aligned boxes)
     * It is assumed that the std::vectors are all of the same length
     * @param rConductivityAreas conductivityAreas[0] is the first region
     * @param rIntraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param rExtraConductivities  extraConductivities[0] is conductivity vector for the first region
     */
    void SetConductivityHeterogeneities(std::vector<ChasteCuboid<3> >& rConductivityAreas,
                                        std::vector< c_vector<double,3> >& rIntraConductivities,
                                        std::vector< c_vector<double,3> >& rExtraConductivities);
    /**
     * Set a number of heterogeneous regions (Axis-aligned ellipsoids)
     * It is assumed that the std::vectors are all of the same length
     * @param conductivityAreas conductivityAreas[0] is the first region
     * @param intraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param extraConductivities  extraConductivities[0] is conductivity vector for the first region
     */
    void SetConductivityHeterogeneitiesEllipsoid(std::vector<ChasteEllipsoid<3> >& conductivityAreas,
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
     *
     * Warning: when specifying output variables, you cannot convert the HDF5
     * output to Meshalyzer, Cmgui or VTK formats, since the converter will get
     * confused by the presence of extra data.  This method thus also turns off
     * visualizer output if the provided vector is non-empty.
     */
    void SetOutputVariables(const std::vector<std::string>& rOutputVariables);

    /**
     * This method may set the output HDF5 file to be written using the original
     * mesh permutation (in situations where a parallel partition may have
     * permuted the node).  The default is to use the new, not original permutation,
     * i.e.  useOriginal=false
     * 
     * @param useOriginal  whether to use the original permutation 
     */
    void SetOutputUsingOriginalNodeOrdering(bool useOriginal);

    /**
     * Set whether the simulation should be checkpointed or not.
     *
     * @param checkpointSimulation whether to do checkpointing
     * @param checkpointTimestep checkpointing timestep
     * @param maxCheckpointsOnDisk maximum number of checkpoint archives to keep on disk
     */
     void SetCheckpointSimulation(bool checkpointSimulation, double checkpointTimestep=-1.0, unsigned maxCheckpointsOnDisk=UINT_MAX);


    // Physiological
    /**
     * 3D version
     * @param rIntraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 3>& rIntraConductivities);
    /**
     * 2D version
     * @param rIntraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 2>& rIntraConductivities);
    /**
     * 1D version
     * @param rIntraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 1>& rIntraConductivities);

    /**
     * 3D version
     * @param rExtraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 3>& rExtraConductivities);
    /**
     * 2D version
     * @param rExtraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 2>& rExtraConductivities);
    /**
     * 1D version
     * @param rExtraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 1>& rExtraConductivities);

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
     * @param kspPreconditioner  a string from {"jacobi", "bjacobi", "hypre", "ml", "spai", "blockdiagonal", "ldufactorisation", "none"}
     */
    void SetKSPPreconditioner(const char* kspPreconditioner);

    /**
     * Set the parameters to be used during mesh adaptation.
     *
     * @param targetError  is the target error passed to the adaptivity library
     * @param sigma  is the value of sigma passed to the adaptivity library
     * @param maxEdgeLength  is the maximum edge length permitted in the adapted mesh
     * @param minEdgeLength  is the minimum edge length permitted in the adapted mesh
     * @param gradation  is the value of gradation passed to the adaptivity library
     * @param maxNodes  is the maximum number of nodes permitted in the adapted mesh
     * @param numSweeps  is the number of adaptive sweeps through the mesh performed by the adaptivity library
     */
    void SetAdaptivityParameters(double targetError, double sigma, double maxEdgeLength, double minEdgeLength,
                                 double gradation, unsigned maxNodes, unsigned numSweeps );

    /**
     * Set the target error to be used during mesh adaptation.
     *
     * @param targetError  is the target error passed to the adaptivity library
     */
    void SetTargetErrorForAdaptivity(double targetError);

    /**
     * Set the value of sigma to be used during mesh adaptation.
     *
     * @param sigma  is the value of sigma passed to the adaptivity library
     */
    void SetSigmaForAdaptivity(double sigma);

    /**
     * Set the maximum edge length to be used during mesh adaptation.
     *
     * @param maxEdgeLength  is the maximum edge length permitted in the adapted mesh
     */
    void SetMaxEdgeLengthForAdaptivity(double maxEdgeLength);

    /**
     * Set the minimum edge length to be used during mesh adaptation.
     *
     * @param minEdgeLength  is the minimum edge length permitted in the adapted mesh
     */
    void SetMinEdgeLengthForAdaptivity(double minEdgeLength);

    /**
     * Set the gradation to be used during mesh adaptation.
     *
     * @param gradation  is the gradation passed to the adaptivity library
     */
    void SetGradationForAdaptivity(double gradation);

    /**
     * Set the maximum number of nodes to be used during mesh adaptation.
     *
     * @param maxNodes  is the maximum number of nodes permitted in the adapted mesh
     */
    void SetMaxNodesForAdaptivity(unsigned maxNodes);

    /**
     * Set the number of adaptive sweeps to be used during mesh adaptation.
     *
     * @param numSweeps  is the number of adaptive sweeps through the mesh performed by the adaptivity library
     */
    void SetNumberOfAdaptiveSweeps(unsigned numSweeps);

    /** Set the parameters of the apd map requested
     *
     *  @param rApdMaps  each entry is a request for a map with
     *  - a percentage in the range [1, 100) (ranges are not checked by this method, but during the calculation)
     *  - a threshold (in mV)
     */
    void SetApdMaps(const std::vector<std::pair<double,double> >& rApdMaps);

    /** Set the parameters of the upstroke time map requested
     *
     *  @param rUpstrokeTimeMaps  is the list of thresholds (? ///\todo improve the description of threshold) with respect to which the upstroke time maps are calculated.
     */
    void SetUpstrokeTimeMaps (std::vector<double>& rUpstrokeTimeMaps);

    /** Set the parameters of the maximal upstroke velocity map requested
     *
     *  @param rMaxUpstrokeVelocityMaps is the list of thresholds (? ///\todo improve the description of threshold) with respect to which the upstroke velocity maps are calculated.
     */
    void SetMaxUpstrokeVelocityMaps (std::vector<double>& rMaxUpstrokeVelocityMaps);

    /** Set the parameters of the conduction velocity map requested
     *
     *  @param rConductionVelocityMaps is a list of origin node indices. One map is created for each origin node.
     */
    void SetConductionVelocityMaps (std::vector<unsigned>& rConductionVelocityMaps);


    // Output visualization

    /** Create the OutputVisualizer element if it doesn't exist */
    void EnsureOutputVisualizerExists(void);

    /** Set whether to convert the output from HDF5 to meshalyzer readable format
     *
     * @param useMeshalyzer
     */
    void SetVisualizeWithMeshalyzer(bool useMeshalyzer=true);

    /** Set whether to convert the output from HDF5 to Cmgui readable format
     *
     * @param useCmgui
     */
    void SetVisualizeWithCmgui(bool useCmgui=true);

    /** Set whether to convert the output from HDF5 to Vtk readable format
     *
     * @param useVtk
     */
    void SetVisualizeWithVtk(bool useVtk=true);
    
    /**
     * Setup electrode parameters.
     *
     *  @param groundSecondElectrode Whether to ground the second electrode (see class documentation)
     *  @param index The value i when applying the electrodes to x_i=a and x_i=b (a<b)
     *  @param magnitude Magnitude of the stimulus
     *  @param startTime Switch on time
     *  @param duration Duration of the stimulus.
     */
    void SetElectrodeParameters( bool groundSecondElectrode,
                                 unsigned index, double magnitude, 
                                 double startTime, double duration );

    /**
     * Get electrode parameters.
     *
     *  @param rGroundSecondElectrode Whether to ground the second electrode (see class documentation)
     *  @param rIndex The value i when applying the electrodes to x_i=a and x_i=b (a<b)
     *  @param rMagnitude Magnitude of the stimulus
     *  @param rStartTime Switch on time
     *  @param rDuration Duration of the stimulus.
     */
    void GetElectrodeParameters(bool& rGroundSecondElectrode,
                                unsigned& rIndex, double& rMagnitude, 
                                double& rStartTime, double& rDuration );
    

private:
    // Only to be accessed by the tests
    friend class TestHeartConfig;

    /*Constructor is private, since the class is only accessed by the singleton instance() method*/
    HeartConfig();

    /** Pointer to parameters read from the user's input XML file
     * (override those given by #mpDefaultParameters).
     */
    boost::shared_ptr<cp::chaste_parameters_type> mpUserParameters;
    /** Pointer to parameters read from the default input XML file (to be read before
     * #mpUserParameters, but may be subsequently overridden).
     */
    boost::shared_ptr<cp::chaste_parameters_type> mpDefaultParameters;

    /** The single instance of the class */
    static std::auto_ptr<HeartConfig> mpInstance;

    /**
     * Whether to read the schema location from the XML file (false) or use the schema
     * located at heart/src/io/ChasteParameters.xsd in the Chaste source tree (true).
     */
    bool mUseFixedSchemaLocation;

    /**
     * Fraction of epicardial layer
     */
    double mEpiFraction;

    /**
     * Fraction of endocardial layer
     */
    double mEndoFraction;

    /**
     * Fraction of midmyocardial layer
     */
    double mMidFraction;

    /**
     * Order index in which the midmyocardial heterogeneities are supplied
     */
    unsigned mIndexMid;

    /**
     * Order index in which the epicardial heterogeneities are supplied
     */
    unsigned mIndexEpi;

    /**
     * Order index in which the endocardial heterogeneities are supplied
     */
    unsigned mIndexEndo;

    /**
     * Flag to check whether the user asked for cellular transmural heterogeneities
     */
    bool mUserAskedForCellularTransmuralHeterogeneities;

   /**
     * Flag to check whether the user asked for cellular heterogeneities with cuboids
     */
    bool mUserAskedForCuboidsForCellularHeterogeneities;

    /**
     * DecideLocation is a convenience method used to get the correct parameter value
     * from the defaults/parameters files.  It checks if the first value  is present and (if not)
     * moves onto the second
     *
     * @param params_ptr  Pointer to quantity within the parameters file (checked first, since it will override a default)
     * @param defaults_ptr  Pointer to quantity within the defaults file (used if there was no override)
     * @param rNameParameter Name of quatity within params_ptr/defaults_ptr (so we can throw a meaningful exception if it's not found)
     */
    template<class TYPE>
    TYPE* DecideLocation(TYPE* params_ptr, TYPE* defaults_ptr, const std::string& rNameParameter) const;

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
     * CheckSimulationIsDefined is a convience method for checking if the "<"ResumeSimulation">" element
     * has been defined and therefore is safe to use the ResumeSimulation().get() pointer to access
     * other data.
     *
     * Throws and exception if not.
     *
     * @param callingMethod string describing the get method performing the check.
     */
     void CheckResumeSimulationIsDefined(std::string callingMethod="") const;

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(HeartConfig)

#endif /*HEARTCONFIG_HPP_*/
