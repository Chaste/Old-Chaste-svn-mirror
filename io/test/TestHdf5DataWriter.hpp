#ifndef TESTHDF5DATAWRITER_HPP_
#define TESTHDF5DATAWRITER_HPP_

#include <hdf5.h>
#include <cxxtest/TestSuite.h>

#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
//#include "Hdf5DataWriter.hpp"

class TestHdf5DataWriter : public CxxTest::TestSuite
{
public:
    void TestSimpleParallelWriteDirectlyWithHdf5()
    {
        // File to write
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory();
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
        hsize_t count[2] = {X, Y};
        hid_t memspace = H5Screate_simple(DIMS, count, NULL);
        
        // Select hyperslab in the file.
        hsize_t offset[2] = {PetscTools::GetMyRank()*X, 0};
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
        
        // Create property list for collective dataset write, and write!  Finally.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
        herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
        TS_ASSERT_EQUALS(status, 0);
        
        // Release resources and close the file
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
        H5Fclose(file_id);
    }
};
#endif /*TESTHDF5DATAWRITER_HPP_*/
