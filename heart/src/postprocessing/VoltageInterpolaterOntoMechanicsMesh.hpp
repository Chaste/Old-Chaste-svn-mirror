#ifndef VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_
#define VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_


#include <vector>
#include <string>
#include "UblasIncludes.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "FineCoarseMeshPair.hpp"
#include "ReplicatableVector.hpp"
#include "HeartConfig.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscTools.hpp"



/**
 *  Very simple one-method class which can be used to convert the voltage from an electrics
 *  (or electromechanics) simulation onto a coarser mechanics mesh, by interpolation. The 
 *  class outputs a HDF5 file corresponding to nodes on the mechanics mesh, and converts it to 
 *  CMGUI output. 
 */ 
template<unsigned DIM>
class VoltageInterpolaterOntoMechanicsMesh
{
public:
    /** 
     *  Constructor, also the main method of the class
     *  
     *  @param rElectricsMesh The electrics mesh
     *  @param rMechanicsMesh The mechanics mesh
     *  @param directory Directory the voltage file is in
     *  @param inputFileNamePrefix Filename (without ".h5") of the electrics solution HDF5 file
     */
    VoltageInterpolaterOntoMechanicsMesh(TetrahedralMesh<DIM,DIM>& rElectricsMesh,
                                         QuadraticMesh<DIM>& rMechanicsMesh,
                                         std::string directory,
                                         std::string inputFileNamePrefix)
    {
        // Read the data from the HDF5 file
        Hdf5DataReader reader(directory,inputFileNamePrefix);

        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
        
        // set up the elements and weights for the coarse nodes in the fine mesh
        FineCoarseMeshPair<DIM> mesh_pair(rElectricsMesh, rMechanicsMesh);
        mesh_pair.SetUpBoxesOnFineMesh();
        mesh_pair.ComputeFineElementsAndWeightsForCoarseNodes(true);
        assert(mesh_pair.rGetElementsAndWeights().size()==rMechanicsMesh.GetNumNodes());

        // create and setup a writer
        Hdf5DataWriter* p_writer = new Hdf5DataWriter(*rMechanicsMesh.GetDistributedVectorFactory(),
                                                      directory,
                                                      "voltage_mechanics_mesh",
                                                      false, //don't clean
                                                      false);

        p_writer->DefineFixedDimension( rMechanicsMesh.GetNumNodes() );
        int voltage_column_id = p_writer->DefineVariable("V","mV");
        p_writer->DefineUnlimitedDimension("Time","msecs");
        p_writer->EndDefineMode();

        // set up a vector to read into
        DistributedVectorFactory factory(rElectricsMesh.GetNumNodes());
        Vec voltage = factory.CreateVec(); 
        std::vector<double> interpolated_voltages(rMechanicsMesh.GetNumNodes());
        Vec voltage_coarse = NULL;

        for(unsigned time_step=0; time_step<num_timesteps; time_step++)
        {
            // read
            reader.GetVariableOverNodes(voltage, "V", time_step);
            ReplicatableVector voltage_repl(voltage);
            
            // interpolate
            for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
            {
                double interpolated_voltage = 0;

                Element<DIM,DIM>& element = *(rElectricsMesh.GetElement(mesh_pair.rGetElementsAndWeights()[i].ElementNum));
                for(unsigned node_index = 0; node_index<element.GetNumNodes(); node_index++)
                {
                    unsigned global_node_index = element.GetNodeGlobalIndex(node_index);
                    interpolated_voltage += voltage_repl[global_node_index]*mesh_pair.rGetElementsAndWeights()[i].Weights(node_index);
                }

                interpolated_voltages[i] = interpolated_voltage;
            }
            
            if(voltage_coarse!=NULL)
            {
                VecDestroy(voltage_coarse);
            }
            voltage_coarse = PetscTools::CreateVec(interpolated_voltages);
          
            // write
            p_writer->PutUnlimitedVariable(time_step);
            p_writer->PutVector(voltage_column_id, voltage_coarse);
            p_writer->AdvanceAlongUnlimitedDimension();            
        }

        if(voltage_coarse!=NULL)
        {
            VecDestroy(voltage);
            VecDestroy(voltage_coarse);
        }

        // delete to flush 
        delete p_writer;
        
        // Convert the new data to CMGUI format.
        // alter the directory in HeartConfig as that is where Hdf5ToCmguiConverter decides
        // where to output
        std::string config_directory = HeartConfig::Instance()->GetOutputDirectory();
        HeartConfig::Instance()->SetOutputDirectory(directory);
        Hdf5ToCmguiConverter<DIM,DIM> converter(directory, 
                                                "voltage_mechanics_mesh", 
                                                &rMechanicsMesh, 
                                                false);
        HeartConfig::Instance()->SetOutputDirectory(config_directory);                                                
    }
};

#endif /*VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_*/
