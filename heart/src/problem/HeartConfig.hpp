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

#include "UblasCustomFunctions.hpp"
#include "ChasteParameters.hpp"
#include <iostream>
#include "Exception.hpp"
#include <vector>
#include "AbstractStimulusFunction.hpp"
#include "SimpleStimulus.hpp"
#include "ChasteCuboid.hpp"
#include "OutputFileHandler.hpp"

#include <boost/shared_ptr.hpp>

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
    void CheckTimeSteps() const;

public:
    /**
     * Call this method to access the global parameters holder.
     *
     * @return a single instance of the class
     */
    static HeartConfig* Instance();

    /**
     * @param fileName The name of the default file - set by default to "ChasteDefaults.xml" on construction
     */
    void SetDefaultsFile(std::string fileName);
    /**
     * #mpUserParameters  is set to a new context associated with a parameters file 
     * @param fileName The name of the parameters file
     */
    void SetParametersFile(std::string fileName);
    /**
     * Write out the complete configuration set as an XML file 
     * Note that the location of ChasteParameters.xsd (schema definition)
     * will be hard-coded in the XML file.
     * @param dirName  Directory path
     * @param fileName Name of file (ideally will have ".xml" suffix)
     */
    void Write(std::string dirName, std::string fileName);
    static void Reset();

    /*
     *  Get methods
     */
    // Simulation
    unsigned GetSpaceDimension() const;
    double GetSimulationDuration() const;
    domain_type GetDomain() const;
    ionic_models_available_type GetDefaultIonicModel() const;
    void GetIonicModelRegions(std::vector<ChasteCuboid>& definedRegions,
                              std::vector<ionic_models_available_type>& ionicModels) const;


    bool GetIsMeshProvided() const;
    bool GetCreateMesh() const;
    bool GetCreateSlab() const;
    bool GetCreateSheet() const;
    bool GetCreateFibre() const;
    bool GetLoadMesh() const;

    void GetSlabDimensions(c_vector<double, 3>& slabDimensions) const;
    void GetSheetDimensions(c_vector<double, 2>& sheetDimensions) const;
    void GetFibreLength(c_vector<double, 1>& fibreLength) const;
    double GetInterNodeSpace() const;

    std::string GetMeshName() const;
    media_type GetConductivityMedia() const;

    void GetStimuli(std::vector<boost::shared_ptr<SimpleStimulus> >& rStimuliApplied, std::vector<ChasteCuboid>& rStimulatedAreas) const;
    void GetCellHeterogeneities(std::vector<ChasteCuboid>& cellHeterogeneityAreas,
                                std::vector<double>& scaleFactorGks,
                                std::vector<double>& scaleFactorIto,
                                std::vector<double>& scaleFactorGkr) const;
    bool GetConductivityHeterogeneitiesProvided() const;
    void GetConductivityHeterogeneities(std::vector<ChasteCuboid>& conductivitiesHeterogeneityAreas,
                                        std::vector< c_vector<double,3> >& intraConductivities,
                                        std::vector< c_vector<double,3> >& extraConductivities) const;
    std::string GetOutputDirectory() const;
    std::string GetOutputFilenamePrefix() const;

    // Physiological
    void GetIntracellularConductivities(c_vector<double, 3>& intraConductivities) const;
    void GetIntracellularConductivities(c_vector<double, 2>& intraConductivities) const;
    void GetIntracellularConductivities(c_vector<double, 1>& intraConductivities) const;

    void GetExtracellularConductivities(c_vector<double, 3>& extraConductivities) const;
    void GetExtracellularConductivities(c_vector<double, 2>& extraConductivities) const;
    void GetExtracellularConductivities(c_vector<double, 1>& extraConductivities) const;

    double GetBathConductivity() const;

    double GetSurfaceAreaToVolumeRatio() const;
    double GetCapacitance() const;

    // Numerical
    double GetOdeTimeStep() const;
    double GetPdeTimeStep() const;
    double GetPrintingTimeStep() const;

    bool GetUseAbsoluteTolerance() const;
    double GetAbsoluteTolerance() const;

    bool GetUseRelativeTolerance() const;
    double GetRelativeTolerance() const;

    const char* GetKSPSolver() const;
    const char* GetKSPPreconditioner() const;

    // Post processing
    bool GetIsPostProcessingRequested() const;
    
    bool GetApdMapsRequested() const;
    void GetApdMaps(std::vector<std::pair<double,double> >& apd_maps) const;
    
    bool GetUpstrokeTimeMapsRequested() const;
    void GetUpstrokeTimeMaps (std::vector<double>& upstroke_time_maps) const;
    
    bool GetIsMaxUpstrokeVelocityMapRequested() const;
    
    bool GetConductionVelocityMapsRequested() const;
    void GetConductionVelocityMaps(std::vector<unsigned>& conduction_velocity_maps) const;


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
     * @param domain  type of simulation bi- mono-domain
     */
    void SetDomain(domain_type domain);
    /**
     * Set the configuration place given cardiac cell models at all mesh nodes 
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
     * @param fibreDefinition  if set (Orthotropic/Axisymmetric) then a fibre file should also be read
     * \todo Is fibre file reading implemented?  If so, where?
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
     * @param outputDirectory  Full path to output directory (will be created if necessary)
     */
    void SetOutputDirectory(std::string outputDirectory);
    /**
     * @param outputFilenamePrefix  Prefix for files
     * If set to "res" this will produce
     * <path>/res.h5
     * <path>/output/res_mesh.pts
     * <path>/output/res_mesh.tri  
     * <path>/output/res_parameters.xml  (a copy of this configuration at the end of the simulation)
     * <path>/output/res_times.info
     * <path>/output/res_V.dat
     */
    void SetOutputFilenamePrefix(std::string outputFilenamePrefix);

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
     * @param kspPreconditioner  a string from {""ilu", "jacobi", "bjacobi", "hypre", "none"}
     */
    void SetKSPPreconditioner(const char* kspPreconditioner);

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
    chaste_parameters_type* mpUserParameters;
    /** Pointer to parameters read from the default input XML file (to be read before 
     * #mpUserParameters, but may be subsequently overridden).
     */
    chaste_parameters_type* mpDefaultParameters;

    /** The single instance of the class */
    static std::auto_ptr<HeartConfig> mpInstance;

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

    /** Utility method to parse an XML parameters file
     * get the parameters using the method 'ChasteParameters(filename)',
     * which returns a std::auto_ptr.
     * @param fileName  Name of XML file
     */    
    chaste_parameters_type* ReadFile(std::string fileName);

};

#endif /*HEARTCONFIG_HPP_*/
