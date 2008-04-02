#ifndef HDF5TOMESHALYZERCONVERTER_HPP_
#define HDF5TOMESHALYZERCONVERTER_HPP_

#include "Hdf5DataReader.hpp"
#include "PetscTools.hpp"

/** 
 *  This class converts from Hdf5 format to meshalyzer format, ie, for
 *  voltage, one file, which looks like
 * 
 *  V_node_0_time_0
 *  ..
 *  V_node_N_time_0
 *  V_node_0_time_1
 *  ..
 *  V_node_N_time_1
 *  V_node_0_time_2
 *  ..
 *  V_node_N_time_M
 * 
 *  The files that are written are <base_name>_V.dat or <base_name>_Phi_e.dat,
 *  where <base_name> is the base name of the original .h5 file. The new files
 *  are written in the same directory as the .h5 file.
 */
class Hdf5ToMeshalyzerConverter
{
private:
    Hdf5DataReader* mpReader;
    std::string mOutputDirectory; 
    std::string mFileBaseName;

    /** A helper method which takes in a string, which must be 'V' or 'Phi_e'
     *  and reads the data corresponding to that string, writing it out in
     *  meshalyzer format.
     */
    void Write(std::string type)
    {
        assert(type=="V" || type=="Phi_e");
        
        OutputFileHandler output_file_handler(mOutputDirectory, false);
        out_stream p_file = output_file_handler.OpenOutputFile(mFileBaseName + "_" + type + ".dat");

        unsigned num_nodes = mpReader->GetNumberOfRows();
        unsigned num_timesteps = mpReader->GetUnlimitedDimensionValues().size();
        
        DistributedVector::SetProblemSize(num_nodes);
        Vec petsc_data_V = DistributedVector::CreateVec();
        for (unsigned time_step=0; time_step<num_timesteps; time_step++)
        {
            mpReader->GetVariableOverNodes(petsc_data_V, type, time_step);
            ReplicatableVector repl_data(petsc_data_V);
            
            assert(repl_data.size()==num_nodes);
            
            if(PetscTools::AmMaster())
            {
                for(unsigned i=0; i<num_nodes; i++)
                {
                    *p_file << repl_data[i] << "\n";
                }
            }
        }
    }


public:
    Hdf5ToMeshalyzerConverter(std::string outputDirectory, 
                              std::string fileBaseName)
    {
        // store dir and filenames, and create a reader
        mOutputDirectory = outputDirectory;
        mFileBaseName = fileBaseName;
        mpReader = new Hdf5DataReader(mOutputDirectory, mFileBaseName);
        
        // check the data file read has one or two variables (ie V; or V and PhiE)
        std::vector<std::string> variable_names = mpReader->GetVariableNames();
        if((variable_names.size()==0) || (variable_names.size()>2))
        {
            EXCEPTION("Data has zero or more than two variables - doesn't appear to be mono or bidomain"); 
        }
        
        // if one variable, a monodomain problem
        if(variable_names.size()==1)
        {
            if(variable_names[0]!="V")
            {
                EXCEPTION("One variable, but it is not called 'V'");
            }
            
            Write("V");
        }

        // if two variable, a bidomain problem
        if(variable_names.size()==2)
        {
            if(variable_names[0]!="V" || variable_names[1]!="Phi_e")
            {
                EXCEPTION("Two variables, but they are not called 'V' and 'Phi_e'");
            }

            Write("V");
            Write("Phi_e");
        }
    }
    
    ~Hdf5ToMeshalyzerConverter()
    {
        delete mpReader;
    }
};

#endif /*HDF5TOMESHALYZERCONVERTER_HPP_*/
