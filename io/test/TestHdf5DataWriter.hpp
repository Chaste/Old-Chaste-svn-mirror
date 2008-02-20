#ifndef TESTHDF5DATAWRITER_HPP_
#define TESTHDF5DATAWRITER_HPP_

#include <hdf5.h>
#include <cxxtest/TestSuite.h>

#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "HDF5DataWriter.hpp"

class TestHdf5DataWriter : public CxxTest::TestSuite
{
public:
    void TestSimpleParallelWriteDirectlyWithHdf5()
    {
        // File to write
        OutputFileHandler oh("hdf5");
        std::string results_dir = oh.GetOutputDirectoryFullPath();
        std::string file_name = results_dir + "test.h5";
        
        // Set up a property list saying how we'll open the file
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL);
        
        // Create a file (collectively) and free the property list
        hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
        
        // Initialise the data for this process
        const unsigned DIMS = 2;
        const unsigned X = 2;
        const unsigned Y = 5;
        int data[X][Y];
        for (unsigned i=0; i<X; i++)
        {
            for (unsigned j=0; j<Y; j++)
            {
                data[i][j] = 100*PetscTools::GetMyRank() + 10*i + j;
            }
        }
        
        // Create the dataspace for the dataset.
        hsize_t dimsf[DIMS]; // dataset dimensions
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        dimsf[0] = X * num_procs;
        dimsf[1] = Y;
        hid_t filespace = H5Screate_simple(DIMS, dimsf, NULL);
        
        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate(file_id, "IntArray", H5T_NATIVE_INT, filespace, H5P_DEFAULT);
        H5Sclose(filespace);

        // Define a dataset in memory for this process
        hsize_t count[DIMS] = {X, Y};
        hid_t memspace = H5Screate_simple(DIMS, count, NULL);
        
        // Select hyperslab in the file.
        hsize_t offset[DIMS] = {PetscTools::GetMyRank()*X, 0};
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
        
        // Create property list for collective dataset write, and write!  Finally.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
        herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
        TS_ASSERT_EQUALS(status, 0);
        
        // Create dataspace for the name, unit attribute
        hsize_t columns[2] = {Y, 21};
        hid_t colspace = H5Screate_simple(1, columns, NULL);
        
        //Create attribute
        char col_data[5][21];
        strcpy(col_data[0], "Noughth");
        strcpy(col_data[1], "First");
        strcpy(col_data[2], "Second");
        strcpy(col_data[3], "Third");
        strcpy(col_data[4], "Fourth");
        
        // create the type 'char'
        hid_t char_type = H5Tcopy(H5T_C_S1);
        //H5Tset_strpad(char_type, H5T_STR_NULLPAD);
        H5Tset_size(char_type, 21 );
        hid_t attr = H5Acreate(dset_id, "Name", char_type, colspace, H5P_DEFAULT  );
        // Write to the attribute        
        status = H5Awrite(attr, char_type, col_data); 
               

        
        //Close dataspace & attribute
        H5Sclose(colspace);
        H5Aclose(attr);
        
        // Release resources and close the file
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
        H5Fclose(file_id);
    }
    
    static const unsigned data_size=17;
    void TestPetscWriteDirectlyWithHdf5()
    {
        
        //Initialise a PETSc vector
        Vec a_vec=PetscTools::CreateVec(data_size);
        double* p_a_vec;
        VecGetArray(a_vec, &p_a_vec);
        int lo, hi;
        VecGetOwnershipRange(a_vec, &lo, &hi);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_a_vec[local_index] = global_index + 100*PetscTools::GetMyRank();
        }
        VecRestoreArray(a_vec, &p_a_vec);
        VecAssemblyBegin(a_vec);
        VecAssemblyEnd(a_vec);
        
        //VecView(a_vec, PETSC_VIEWER_STDOUT_WORLD);

        // File to write
        OutputFileHandler oh("hdf5", false);
        std::string results_dir = oh.GetOutputDirectoryFullPath();
        std::string file_name = results_dir + "vec.h5";
        
        // Set up a property list saying how we'll open the file
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL);
        
        // Create a file (collectively) and free the property list
        hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
        
        const unsigned DIMS = 1;
        
        // Create the dataspace for the dataset.
        //TS_ASSERT_EQUALS(data_size, hi-lo);
        hsize_t dimsf[DIMS]={data_size}; // dataset dimensions
        
        hid_t filespace = H5Screate_simple(DIMS, dimsf, NULL);
        
        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate(file_id, "TheVector", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT);
        H5Sclose(filespace);

        // Define a dataset in memory for this process
        hsize_t count[DIMS] = {hi-lo};
        hid_t memspace = H5Screate_simple(DIMS, count, NULL);
        
        // Select hyperslab in the file.
        hsize_t offset[DIMS] = {lo};
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
        
        // Create property list for collective dataset write, and write!  Finally.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
        VecGetArray(a_vec, &p_a_vec);
        herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, p_a_vec);
        VecRestoreArray(a_vec, &p_a_vec);
        
        TS_ASSERT_EQUALS(status, 0);
        
        // Release resources and close the file
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
        H5Fclose(file_id);
        
        VecDestroy(a_vec);
    }
    
    void TestReadAndChecksumHdf5()
    { 
        // File to read
        OutputFileHandler oh("hdf5", false);
        double data[data_size];
        std::string results_dir = oh.GetOutputDirectoryFullPath();
        std::string file_name = results_dir + "vec.h5";
        
        hsize_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        //H5Tget_nmembers(file_id);
        hsize_t dataset_id = H5Dopen(file_id, "TheVector");
        hsize_t dxpl = H5Pcreate(H5P_DATASET_XFER);
        hsize_t edc = H5Pget_edc_check(dxpl);
        TS_ASSERT_EQUALS(edc, (hsize_t) 1) //Checksum is enabled
        
        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl, 
            data);
       
        TS_ASSERT_EQUALS(status, 0); 

        //Check the index
        for (unsigned i=0;i<data_size;i++)
        {
            TS_ASSERT_EQUALS(((unsigned)data[i]%100), i);
        }
        //Check the final component
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        //The last component was owned by processor "num_procs-1"
        TS_ASSERT_EQUALS(((int)data[data_size-1]/100), num_procs-1);
        
        H5Pclose (dxpl);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
    }
  
    void TestHDF5DataWriterPut() throw(Exception)
    {
        HDF5DataWriter writer("", "hdf5_test");
        writer.DefineFixedDimension("Node","dimensionless",4);
       // int ik_id = 
        writer.DefineVariable("I_K","milliamperes");

//// get an error in the H5Pset_fapl_mpio() call in the following method
//        writer.EndDefineMode();
//
//        writer.PutVariable(ik_id,  0.0, 0);
//        writer.PutVariable(ik_id, 10.0, 1);
//        writer.PutVariable(ik_id, 20.0, 2);
//        writer.PutVariable(ik_id, 30.0, 3);
        
        writer.Close();
        
    }
};
#endif /*TESTHDF5DATAWRITER_HPP_*/
