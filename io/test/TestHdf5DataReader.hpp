#ifndef TESTHDF5READER_HPP_
#define TESTHDF5READER_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "Hdf5DataWriter.hpp"
#include "Hdf5DataReader.hpp"

class TestHdf5DataReader : public CxxTest::TestSuite
{
private :
    #define H5FILE_NAME "SDS.h5"
    #define DATASETNAME "IntArray"

    void WriteDataTestSimpleReadDirectlyWithHdf5()
    {
        int const NX = 5;                      /* dataset dimensions */
        int const NY = 6;
        int const RANK = 2;

        hid_t       file, dataset;         /* file and dataset handles */
        hid_t       datatype, dataspace;   /* handles */
        hsize_t     dimsf[2];              /* dataset dimensions */
        herr_t      status;
        int         data[NX][NY];          /* data to write */
        int         i, j;
    
        /*
         * Data  and output buffer initialization.
         */
        for (j = 0; j < NX; j++) {
        for (i = 0; i < NY; i++)
            data[j][i] = i + j;
        }
        /*
         * 0 1 2 3 4 5
         * 1 2 3 4 5 6
         * 2 3 4 5 6 7
         * 3 4 5 6 7 8
         * 4 5 6 7 8 9
         */
    
        /*
         * Create a new file using H5F_ACC_TRUNC access,
         * default file creation properties, and default file
         * access properties.
         */
        file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
        /*
         * Describe the size of the array and create the data space for fixed
         * size dataset.
         */
        dimsf[0] = NX;
        dimsf[1] = NY;
        dataspace = H5Screate_simple(RANK, dimsf, NULL);
    
        /*
         * Define datatype for the data in the file.
         * We will store little endian INT numbers.
         */
        datatype = H5Tcopy(H5T_NATIVE_INT);
        status = H5Tset_order(datatype, H5T_ORDER_LE);
    
        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace,
                H5P_DEFAULT);
    
        /*
         * Write the data to the dataset using default transfer properties.
         */
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, data);
    
        /*
         * Close/release resources.
         */
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Fclose(file);
        
    }

public :

    void TestSimpleReadDirectlyWithHdf5()
    {
        
        int const NX_SUB = 3;           /* hyperslab dimensions */ 
        int const NY_SUB = 4; 
        int const NX = 7;           /* output buffer dimensions */ 
        int const NY = 7; 
        int const NZ = 3; 
        //int const RANK      =   2;
        int const RANK_OUT  =   3;
        
        hid_t       file, dataset;         /* handles */
        hid_t       datatype, dataspace;   
        hid_t       memspace; 
        H5T_class_t datatype_class;                 /* datatype class */
        H5T_order_t order;                 /* data order */
        size_t      size;                  /*
                            * size of the data element         
                            * stored in file
                            */
        hsize_t     dimsm[3];              /* memory space dimensions */
        hsize_t     dims_out[2];           /* dataset dimensions */      
        herr_t      status;                             
    
        int         data_out[NX][NY][NZ ]; /* output buffer */
       
        hsize_t      count[2];              /* size of the hyperslab in the file */
        hsize_t      offset[2];             /* hyperslab offset in the file */
        hsize_t      count_out[3];          /* size of the hyperslab in memory */
        hsize_t      offset_out[3];         /* hyperslab offset in memory */
        int          i, j, k, status_n, rank;
    
    
        // Create the file it's gonna be read
        WriteDataTestSimpleReadDirectlyWithHdf5();
        
        for (j = 0; j < NX; j++) 
        {
            for (i = 0; i < NY; i++) 
            {
                for (k = 0; k < NZ ; k++)
                    data_out[j][i][k] = 0;
            }
        } 
     
        /*
         * Open the file and the dataset.
         */
        file = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
        dataset = H5Dopen(file, DATASETNAME);
    
        /*
         * Get datatype and dataspace handles and then query
         * dataset class, order, size, rank and dimensions.
         */
        datatype  = H5Dget_type(dataset);     /* datatype handle */ 
        datatype_class     = H5Tget_class(datatype);
        //if (datatype_class == H5T_INTEGER) printf("Data set has INTEGER type \n");
        order     = H5Tget_order(datatype);
        //if (order == H5T_ORDER_LE) printf("Little endian order \n");
    
        size  = H5Tget_size(datatype);
        //printf(" Data size is %d \n", size);
    
        dataspace = H5Dget_space(dataset);    /* dataspace handle */
        rank      = H5Sget_simple_extent_ndims(dataspace);
        status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
        //printf("rank %d, dimensions %lu x %lu \n", rank,
        //   (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));
    
        /* 
         * Define hyperslab in the dataset. 
         */
        offset[0] = 1;
        offset[1] = 2;
        count[0]  = NX_SUB;
        count[1]  = NY_SUB;
        status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, 
                     count, NULL);
    
        /*
         * Define the memory dataspace.
         */
        dimsm[0] = NX;
        dimsm[1] = NY;
        dimsm[2] = NZ ;
        memspace = H5Screate_simple(RANK_OUT,dimsm,NULL);   
    
        /* 
         * Define memory hyperslab. 
         */
        offset_out[0] = 3;
        offset_out[1] = 0;
        offset_out[2] = 0;
        count_out[0]  = NX_SUB;
        count_out[1]  = NY_SUB;
        count_out[2]  = 1;
        status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, 
                     count_out, NULL);
    
        /*
         * Read data from hyperslab in the file into the hyperslab in 
         * memory and display.
         */
        status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace,
                 H5P_DEFAULT, data_out);

        /*
         * 0 0 0 0 0 0 0
         * 0 0 0 0 0 0 0
         * 0 0 0 0 0 0 0
         * 3 4 5 6 0 0 0  
         * 4 5 6 7 0 0 0
         * 5 6 7 8 0 0 0
         * 0 0 0 0 0 0 0
         */
        for (int row=0; row<NX; row++)
        {
            for (int column=0; column<NY; column++)
            {
                if (row>2 && row<6 && column<4) 
                {
                    TS_ASSERT_EQUALS(data_out[row][column][0], row+column);
                }
                else
                {
                    TS_ASSERT_EQUALS(data_out[row][column][0], 0);
                }    
            }
        }
    
    
        /*
         * Close/release resources.
         */
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Sclose(memspace);
        H5Fclose(file);
    }     

private:
    const static unsigned number_nodes = 100;

    void WriteMultiStepData()
    {       
        DistributedVector::SetProblemSize(number_nodes);
               
        Hdf5DataWriter writer("hdf5_reader", "hdf5_test_complete_format", false);
        writer.DefineFixedDimension(number_nodes);
        
        int node_id = writer.DefineVariable("Node","dimensionless");
        int ik_id = writer.DefineVariable("I_K","milliamperes");
        int ina_id = writer.DefineVariable("I_Na","milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();
      
        Vec petsc_data_1=DistributedVector::CreateVec();
        DistributedVector distributed_vector_1(petsc_data_1);
        
        Vec petsc_data_2=DistributedVector::CreateVec();
        DistributedVector distributed_vector_2(petsc_data_2);

        Vec petsc_data_3=DistributedVector::CreateVec();
        DistributedVector distributed_vector_3(petsc_data_3);
        
        for (unsigned time_step=0; time_step<10; time_step++)
        {
            // write some values
            for (DistributedVector::Iterator index = DistributedVector::Begin();
                 index!= DistributedVector::End();
                 ++index)
            {
                distributed_vector_1[index] =  index.Global;
                distributed_vector_2[index] =  time_step*1000 + 100 + index.Global;
                distributed_vector_3[index] =  time_step*1000 + 200 + index.Global;
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();
                
            // write the vector
            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(ik_id, petsc_data_2);
            writer.PutVector(ina_id, petsc_data_3);
            writer.PutUnlimitedVariable(time_step);
            
            // Don't advance in the last iteration
            if (time_step < 9)
            {
                writer.AdvanceAlongUnlimitedDimension();
            }
        }
        
        writer.Close();
    }


public:
    void TestMultiStepReader() throw (Exception)
    {
        WriteMultiStepData();
        
        Hdf5DataReader reader("hdf5_reader", "hdf5_test_complete_format");

        for (unsigned node_index=0; node_index<number_nodes; node_index++)
        {        
            std::vector<double> node_values = reader.GetVariableOverTime("Node", node_index);
            std::vector<double> i_k_values = reader.GetVariableOverTime("I_K", node_index);
            std::vector<double> i_na_values = reader.GetVariableOverTime("I_Na", node_index);

            TS_ASSERT_EQUALS(node_values.size(), 10u);
            TS_ASSERT_EQUALS(i_k_values.size(), 10u);
            TS_ASSERT_EQUALS(i_na_values.size(), 10u);

            for(unsigned i=0; i<node_values.size(); i++)
            {
                TS_ASSERT_DELTA( node_values[i], node_index, 1e-9);
                TS_ASSERT_DELTA( i_k_values[i], i*1000 + 100 + node_index, 1e-9);
                TS_ASSERT_DELTA( i_na_values[i], i*1000 + 200 + node_index, 1e-9);               
            }
        }

        int number_nodes=100;
        DistributedVector::SetProblemSize(number_nodes);

        Vec petsc_data_1=DistributedVector::CreateVec();
        DistributedVector distributed_vector_1(petsc_data_1);
        
        Vec petsc_data_2=DistributedVector::CreateVec();
        DistributedVector distributed_vector_2(petsc_data_2);

        Vec petsc_data_3=DistributedVector::CreateVec();
        DistributedVector distributed_vector_3(petsc_data_3);
        
        for (unsigned time_step=0; time_step<10; time_step++)
        {
            reader.GetVariableOverNodes(petsc_data_1, "Node", time_step);
            reader.GetVariableOverNodes(petsc_data_2, "I_K", time_step);
            reader.GetVariableOverNodes(petsc_data_3, "I_Na", time_step);
            
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();
            
 
            // check values
            for (DistributedVector::Iterator index = DistributedVector::Begin();
                 index!= DistributedVector::End();
                 ++index)
            {
                TS_ASSERT_EQUALS(distributed_vector_1[index], index.Global);
                TS_ASSERT_EQUALS(distributed_vector_2[index], time_step*1000 + 100 + index.Global);
                TS_ASSERT_EQUALS(distributed_vector_3[index], time_step*1000 + 200 + index.Global);
            }
                
        }

        std::vector<double> unlimited_values = reader.GetUnlimitedDimensionValues();
        
        for (unsigned i=0; i< unlimited_values.size(); i++)
        {
            TS_ASSERT_EQUALS(unlimited_values[i], i); 
        }

        reader.Close();
    }
    
    void TestNonMultiStepExceptions ()
    {
        DistributedVector::SetProblemSize(number_nodes);
               
        Hdf5DataWriter writer("hdf5_reader", "hdf5_test_overtime_exceptions", false);
        writer.DefineFixedDimension(number_nodes);
        
        writer.DefineVariable("Node","dimensionless");
        writer.DefineVariable("I_K","milliamperes");
        writer.DefineVariable("I_Na","milliamperes");

        writer.EndDefineMode();
        writer.Close();
     
        Hdf5DataReader reader("hdf5_reader", "hdf5_test_overtime_exceptions");
        
        TS_ASSERT_THROWS_ANYTHING(reader.GetVariableOverTime("Node", 99/*node*/));
       
        Vec data=DistributedVector::CreateVec();
        TS_ASSERT_THROWS_ANYTHING(reader.GetVariableOverNodes(data, "Node", 1/*timestep*/));        
                
        reader.Close();        
    }

    void TestMultiStepExceptions ()
    {
        DistributedVector::SetProblemSize(number_nodes);
               
        Hdf5DataWriter writer("hdf5_reader", "hdf5_test_overtime_exceptions", false);
        writer.DefineFixedDimension(number_nodes);
        
        writer.DefineVariable("Node","dimensionless");
        writer.DefineVariable("I_K","milliamperes");
        writer.DefineVariable("I_Na","milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();
        
        writer.AdvanceAlongUnlimitedDimension();
        
        writer.Close();
     
        TS_ASSERT_THROWS_ANYTHING(Hdf5DataReader reader2("hdf5_reader", "hdf5_wrong_name"));     
        Hdf5DataReader reader("hdf5_reader", "hdf5_test_overtime_exceptions");
               
        TS_ASSERT_THROWS_NOTHING(reader.GetVariableOverTime("Node", 99/*node*/));
        TS_ASSERT_THROWS_ANYTHING(reader.GetVariableOverTime("WrongName", 99/*node*/));
        TS_ASSERT_THROWS_ANYTHING(reader.GetVariableOverTime("Node", 100/*node*/));
               
        Vec data=DistributedVector::CreateVec();
        TS_ASSERT_THROWS_NOTHING(reader.GetVariableOverNodes(data, "Node", 1/*timestep*/));
        TS_ASSERT_THROWS_ANYTHING(reader.GetVariableOverNodes(data, "WrongName"));
        TS_ASSERT_THROWS_ANYTHING(reader.GetVariableOverNodes(data, "I_K", 2/*timestep*/));
                
        reader.Close();        
    }

};

#endif /*TESTHDF5READER_HPP_*/
