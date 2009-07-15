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
#ifndef TESTHDF5DATAWRITER_HPP_
#define TESTHDF5DATAWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "Hdf5DataWriter.hpp"
#include "Hdf5DataReader.hpp"
#include "DistributedVectorFactory.hpp"

class TestHdf5DataWriter : public CxxTest::TestSuite
{
private:

    Hdf5DataWriter* mpTestWriter;

    bool CompareFilesViaHdf5DataReader(std::string pathname1, std::string filename1, bool makeAbsolute1,
        std::string pathname2, std::string filename2, bool makeAbsolute2)
    {
        Hdf5DataReader reader1(pathname1, filename1, makeAbsolute1);
        Hdf5DataReader reader2(pathname2, filename2, makeAbsolute2);

        unsigned number_nodes1 = reader1.GetNumberOfRows();
        unsigned number_nodes2 = reader2.GetNumberOfRows();
        if (number_nodes1 != number_nodes2)
        {
            std::cout << "Number of nodes " << number_nodes1 << " and " << number_nodes2 << " don't match\n";
            return false;
        }
        // Check the variable names and units
        std::vector<std::string> variable_names1 = reader1.GetVariableNames();
        std::vector<std::string> variable_names2 = reader2.GetVariableNames();
        unsigned num_vars = variable_names1.size();
        if (num_vars != variable_names2.size())
        {
            std::cout << "Number of variables " << variable_names1.size()
                      << " and " << variable_names2.size() << " don't match\n";
            return false;
        }
        for (unsigned var=0; var<num_vars; var++)
        {
            std::string var_name = variable_names1[var];
            if (var_name != variable_names2[var])
            {
                std::cout << "Variable names " << var_name << " and "
                          << variable_names2[var] << " don't match\n";
                return false;
            }
            if (reader1.GetUnit(var_name) != reader2.GetUnit(var_name))
            {
                std::cout << "Units names " << reader1.GetUnit(var_name)
                          << " and " << reader2.GetUnit(var_name) << " don't match\n";
                return false;
            }
        }
        // Check the timestep vectors
        std::vector<double> times1 = reader1.GetUnlimitedDimensionValues();
        std::vector<double> times2 = reader2.GetUnlimitedDimensionValues();

        if (times1.size() != times2.size())
        {
            std::cout << "Time step sizes " << times1.size()
                      << " and " << times2.size() << " don't match\n";
            return false;
        }

        for (unsigned timestep=0; timestep<times1.size(); timestep++)
        {
            if (times1[timestep]!=times2[timestep])
            {
                std::cout << "Time steps " << times1[timestep]
                          << " and " << times2[timestep] << " don't match\n";
                return false;
            }
        }

        bool is_complete1 = reader1.IsDataComplete();
        bool is_complete2 = reader2.IsDataComplete();

        if (is_complete1 != is_complete2)
        {
            std::cout<<"One of the readers has incomplete data and the other doesn't\n";
            return false;
        }

        if (is_complete1)
        {
            DistributedVectorFactory factory(number_nodes1);

            Vec data1 = factory.CreateVec();
            Vec data2 = factory.CreateVec();

            for (unsigned timestep=0; timestep<times1.size(); timestep++)
            {
                for (unsigned var=0; var<num_vars; var++)
                {
                    PetscTruth is_equal;
                    reader1.GetVariableOverNodes(data1, variable_names1[var], timestep);
                    reader2.GetVariableOverNodes(data2, variable_names2[var], timestep);
                    VecEqual(data1, data2, &is_equal);
                    /*
                     *
                    std::cout<<"timestep "<<timestep<< "variable_name "<< variable_names1[var]<<"\n";
                    std::cout<<"First------\n";
                    VecView(data1, 0);
                    std::cout<<"Second-----\n";
                    VecView(data1, 0);
                    */

                    if (is_equal != PETSC_TRUE)
                    {
                        return false;
                    }
                }
            }
           VecDestroy(data1);
           VecDestroy(data2);
        }
        else
        {
            // Incomplete data

            // Check the index vectors
            std::vector<unsigned> indices1 = reader1.GetIncompleteNodeMap();
            std::vector<unsigned> indices2 = reader2.GetIncompleteNodeMap();

            if (indices1.size() != indices2.size())
            {
                std::cout << "Index map sizes " << indices1.size() << " and " << indices2.size() << " don't match\n";
                return false;
            }

            for (unsigned index=0; index<indices1.size(); index++)
            {
                if (indices1[index]!=indices2[index])
                {
                   std::cout << "Time steps " << indices1[index] << " and " << indices2[index] << " don't match\n";
                   return false;
                }
            }

            // Check all the data
            for (unsigned index=0; index<indices1.size(); index++)
            {
                unsigned node_index = indices1[index];
                for (unsigned var=0; var<num_vars; var++)
                {
                  std::vector<double> var_over_time1 = reader1.GetVariableOverTime(variable_names1[var], node_index);
                  std::vector<double> var_over_time2 = reader2.GetVariableOverTime(variable_names1[var], node_index);
                  for (unsigned time_step=0;time_step< var_over_time1.size(); time_step++)
                  {
                     if (var_over_time1[time_step] != var_over_time2[time_step])
                     {
                        std::cout<<"Node "<<node_index<<" at time step "<<time_step<<" variable "<<variable_names1[var]<<
                            " differs ("<<var_over_time1[time_step]<<" != "<<var_over_time2[time_step]<<")\n";
                     }
                  }
                }
            }
        }
       return true;
    }

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

        // Create attribute
        char col_data[5][21];
        strcpy(col_data[0], "Noughth");
        strcpy(col_data[1], "First");
        strcpy(col_data[2], "Second");
        strcpy(col_data[3], "Third");
        strcpy(col_data[4], "Fourth");

        // Create the type 'char'
        hid_t char_type = H5Tcopy(H5T_C_S1);
        // H5Tset_strpad(char_type, H5T_STR_NULLPAD);
        H5Tset_size(char_type, 21 );
        hid_t attr = H5Acreate(dset_id, "Name", char_type, colspace, H5P_DEFAULT  );
        // Write to the attribute
        status = H5Awrite(attr, char_type, col_data);

        // Close dataspace & attribute
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
        // Initialise a PETSc vector
        Vec a_vec = PetscTools::CreateVec(data_size);
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

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl, data);

        TS_ASSERT_EQUALS(status, 0);

        // Check the index
        for (unsigned i=0; i<data_size; i++)
        {
            TS_ASSERT_EQUALS(((unsigned)data[i]%100), i);
        }

        // Check the final component
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        // The last component was owned by processor "num_procs-1"
        TS_ASSERT_EQUALS(((int)data[data_size-1]/100), num_procs-1);

        H5Pclose (dxpl);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
    }

    void TestHdf5DataWriterMultipleColumns() throw(Exception)
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "hdf5", "hdf5_test_multi_column", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node","dimensionless");
        int ik_id = writer.DefineVariable("I_K","milliamperes");
        int ina_id = writer.DefineVariable("I_Na","milliamperes");

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        // Write some values
        for (DistributedVector::Iterator index = distributed_vector_1.Begin();
             index!= distributed_vector_1.End();
             ++index)
        {
            distributed_vector_1[index] =  index.Global;
            distributed_vector_2[index] =  1000 + index.Global;
        }
        distributed_vector_1.Restore();
        distributed_vector_2.Restore();

        // Write the vector
        writer.PutVector(node_id, petsc_data_1);
        writer.PutVector(ik_id, petsc_data_1);
        writer.PutVector(ina_id, petsc_data_2);

        writer.Close();

//        if (PetscTools::AmMaster())
//        {
//            // call h5dump to take the binary hdf5 output file and print it
//            // to a text file. Note that the first line of the txt file would
//            // be the directory it has been printed to, but is this line is
//            // removed by piping the output through sed to delete the first line
//            OutputFileHandler handler("hdf5",false);
//            std::string file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_multi_column.h5";
//            std::string new_file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_multi_column_dumped.txt";
//            system( ("h5dump "+file+" | sed 1d > "+new_file).c_str() );
//
//            TS_ASSERT_EQUALS(system(("diff " + new_file + " io/test/data/hdf5_test_multi_column_dumped.txt").c_str()), 0);
//        }

        TS_ASSERT(CompareFilesViaHdf5DataReader("hdf5", "hdf5_test_multi_column", true,
            "io/test/data", "hdf5_test_multi_column", false));

        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);
    }

    void TestHdf5DataWriterNonEvenRowDistribution() throw(Exception)
    {
        int number_nodes = 100;

        PetscInt local_number_of_nodes;

        if (PetscTools::AmMaster())
        {
            local_number_of_nodes = number_nodes - PetscTools::GetNumProcs() + 1;
        }
        else
        {
            local_number_of_nodes = 1;
        }

        DistributedVectorFactory factory(number_nodes, local_number_of_nodes);

        Hdf5DataWriter writer(factory, "hdf5", "hdf5_non_even_row_dist", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node","dimensionless");
        int ik_id = writer.DefineVariable("I_K","milliamperes");
        int ina_id = writer.DefineVariable("I_Na","milliamperes");

        writer.EndDefineMode();

        Vec petsc_data_1=factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2=factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        // Write some values
        for (DistributedVector::Iterator index = distributed_vector_1.Begin();
             index!= distributed_vector_1.End();
             ++index)
        {
            distributed_vector_1[index] = index.Global;
            distributed_vector_2[index] = 1000 + index.Global;
        }
        distributed_vector_1.Restore();
        distributed_vector_2.Restore();

        // Write the vector
        writer.PutVector(node_id, petsc_data_1);
        writer.PutVector(ik_id, petsc_data_1);
        writer.PutVector(ina_id, petsc_data_2);

        writer.Close();

//        if (PetscTools::AmMaster())
//        {
//            // call h5dump to take the binary hdf5 output file and print it
//            // to a text file. Note that the first line of the txt file would
//            // be the directory it has been printed to, but is this line is
//            // removed by piping the output through sed to delete the first line
//            OutputFileHandler handler("hdf5",false);
//            std::string file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_multi_column.h5";
//            std::string new_file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_multi_column_dumped.txt";
//            system( ("h5dump "+file+" | sed 1d > "+new_file).c_str() );
//
//            TS_ASSERT_EQUALS(system(("diff " + new_file + " io/test/data/hdf5_test_multi_column_dumped.txt").c_str()), 0);
//        }

        TS_ASSERT(CompareFilesViaHdf5DataReader("hdf5", "hdf5_non_even_row_dist", true,
                                                "io/test/data", "hdf5_test_multi_column", false));

        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);
    }


    void TestHdf5DataWriterFullFormatIncomplete() throw(Exception)
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "hdf5", "hdf5_test_full_format_incomplete", false);

        int node_id = writer.DefineVariable("Node","dimensionless");
        int ik_id = writer.DefineVariable("I_K","milliamperes");
        int ina_id = writer.DefineVariable("I_Na","milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(60);
        writer.DefineFixedDimension(node_numbers, number_nodes);

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3 = factory.CreateVec();
        DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_3);

        for (unsigned time_step=0; time_step<10; time_step++)
        {
            // Write some values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index!= distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] =  index.Global;
                distributed_vector_2[index] =  time_step*1000 + 100 + index.Global;
                distributed_vector_3[index] =  time_step*1000 + 200 + index.Global;
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();

            // Write the vector

            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(ik_id, petsc_data_2);
            writer.PutVector(ina_id, petsc_data_3);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

//        if (PetscTools::AmMaster())
//        {
//            // call h5dump to take the binary hdf5 output file and print it
//            // to a text file. Note that the first line of the txt file would
//            // be the directory it has been printed to, but is this line is
//            // removed by piping the output through sed to delete the first line
//            OutputFileHandler handler("hdf5",false);
//            std::string file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format_incomplete.h5";
//            std::string new_file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format_incomplete_dumped.txt";
//            system( ("h5dump "+file+" | sed 1d > "+new_file).c_str() );
//
//            TS_ASSERT_EQUALS(system(("diff " + new_file + " io/test/data/hdf5_test_full_format_incomplete_dumped.txt").c_str()), 0);
//        }

        TS_ASSERT(CompareFilesViaHdf5DataReader("hdf5", "hdf5_test_full_format_incomplete", true,
                                                "io/test/data", "hdf5_test_full_format_incomplete", false));

        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);
        VecDestroy(petsc_data_3);
    }

    void TestHdf5DataWriterFullFormat() throw(Exception)
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "hdf5", "hdf5_test_full_format", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int ik_id = writer.DefineVariable("I_K", "milliamperes");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3 = factory.CreateVec();
        DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_3);

        for (unsigned time_step=0; time_step<10; time_step++)
        {
            // Write some values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index!= distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = index.Global;
                distributed_vector_2[index] = time_step*1000 + 100 + index.Global;
                distributed_vector_3[index] = time_step*1000 + 200 + index.Global;
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();

            // Write the vector
            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(ik_id, petsc_data_2);
            writer.PutVector(ina_id, petsc_data_3);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

//        if (PetscTools::AmMaster())
//        {
//            // call h5dump to take the binary hdf5 output file and print it
//            // to a text file. Note that the first line of the txt file would
//            // be the directory it has been printed to, but is this line is
//            // removed by piping the output through sed to delete the first line
//            OutputFileHandler handler("hdf5",false);
//            std::string file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format.h5";
//            std::string new_file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format_dumped.txt";
//            system( ("h5dump "+file+" | sed 1d > "+new_file).c_str() );
//
//            TS_ASSERT_EQUALS(system(("diff " + new_file + " io/test/data/hdf5_test_full_format_dumped.txt").c_str()), 0);
//        }
         TS_ASSERT(CompareFilesViaHdf5DataReader("hdf5", "hdf5_test_full_format", true,
                                                 "io/test/data", "hdf5_test_full_format", false));

        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);
        VecDestroy(petsc_data_3);
    }

    void TestHdf5DataWriterFullFormatStriped() throw(Exception)
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory, "hdf5", "hdf5_test_full_format_striped", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");

        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        DistributedVectorFactory factory(number_nodes);

        Vec petsc_data_short = vec_factory.CreateVec();
        DistributedVector distributed_vector_short = vec_factory.CreateDistributedVector(petsc_data_short);

        Vec node_number = vec_factory.CreateVec();
        DistributedVector distributed_node_number = vec_factory.CreateDistributedVector(node_number);

        for (DistributedVector::Iterator index = distributed_vector_short.Begin();
             index!= distributed_vector_short.End();
             ++index)
        {
            distributed_node_number[index] = index.Global;
            distributed_vector_short[index] = -0.5;
        }
        distributed_node_number.Restore();
        distributed_vector_short.Restore();

        Vec petsc_data_long = factory.CreateVec(2);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector_long,1 );

        for (unsigned time_step=0; time_step<10; time_step++)
        {
            for (DistributedVector::Iterator index = distributed_vector_long.Begin();
                 index!= distributed_vector_long.End();
                 ++index)
            {
                vm_stripe[index] =  time_step*1000 + index.Global*2;
                phi_e_stripe[index] =  time_step*1000 + index.Global*2+1;
            }
            distributed_vector_long.Restore();

            writer.PutVector(node_id, node_number);
            writer.PutVector(ina_id, petsc_data_short);
            writer.PutStripedVector(vm_id, phi_e_id, petsc_data_long);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

//        if (PetscTools::AmMaster())
//        {
//            // call h5dump to take the binary hdf5 output file and print it
//            // to a text file. Note that the first line of the txt file would
//            // be the directory it has been printed to, but is this line is
//            // removed by piping the output through sed to delete the first line
//            OutputFileHandler handler("hdf5",false);
//            std::string file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format_striped.h5";
//            std::string new_file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format_striped_dumped.txt";
//            system( ("h5dump "+file+" | sed 1d > "+new_file).c_str() );
//
//            TS_ASSERT_EQUALS(system(("diff " + new_file + " io/test/data/hdf5_test_full_format_striped_dumped.txt").c_str()), 0);
//        }
//

        TS_ASSERT(CompareFilesViaHdf5DataReader("hdf5", "hdf5_test_full_format_striped", true,
                                                "io/test/data", "hdf5_test_full_format_striped", false));

        VecDestroy(node_number);
        VecDestroy(petsc_data_long);
        VecDestroy(petsc_data_short);
    }

    void TestHdf5DataWriterFullFormatStripedIncomplete() throw(Exception)
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory, "hdf5", "hdf5_test_full_format_striped_incomplete", false);

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(60);
        writer.DefineFixedDimension(node_numbers, number_nodes);

        int vm_id = writer.DefineVariable("V_m","millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e","millivolts");

        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        DistributedVectorFactory factory(number_nodes);
        Vec petsc_data_long = factory.CreateVec(2);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector_long,1 );

        for (unsigned time_step=0; time_step<2; time_step++)
        {
            for (DistributedVector::Iterator index = distributed_vector_long.Begin();
                 index!= distributed_vector_long.End();
                 ++index)
            {
                vm_stripe[index] = (time_step+1)*1000 + index.Global;
                phi_e_stripe[index] =  index.Global;
            }
            distributed_vector_long.Restore();

            writer.PutStripedVector(vm_id, phi_e_id, petsc_data_long);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

//        if (PetscTools::AmMaster())
//        {
//            // call h5dump to take the binary hdf5 output file and print it
//            // to a text file. Note that the first line of the txt file would
//            // be the directory it has been printed to, but is this line is
//            // removed by piping the output through sed to delete the first line
//            OutputFileHandler handler("hdf5",false);
//            std::string file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format_striped_incomplete.h5";
//            std::string new_file = handler.GetOutputDirectoryFullPath() + "/hdf5_test_full_format_striped_incomplete_dumped.txt";
//            system( ("h5dump "+file+" | sed 1d > "+new_file).c_str() );
//
//            TS_ASSERT_EQUALS(system(("diff " + new_file + " io/test/data/hdf5_test_full_format_striped_incomplete_dumped.txt").c_str()), 0);
//        }
//

        TS_ASSERT(CompareFilesViaHdf5DataReader("hdf5", "hdf5_test_full_format_striped_incomplete", true,
                                                "io/test/data", "hdf5_test_full_format_striped_incomplete", false));

        VecDestroy(petsc_data_long);
    }

    void TestNonImplementedFeatures()
    {
        int number_nodes = 100;
        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "hdf5", "hdf5_test_non_implemented", false);
        writer.DefineFixedDimension(number_nodes);

        int vm_id = writer.DefineVariable("V_m","millivolts");
        int ina_id = writer.DefineVariable("I_Na","milliamperes");
        int phi_e_id = writer.DefineVariable("Phi_e","millivolts");

        writer.EndDefineMode();

        Vec petsc_data_short = factory.CreateVec();
        DistributedVector distributed_vector_short = factory.CreateDistributedVector(petsc_data_short);

        for (DistributedVector::Iterator index = distributed_vector_short.Begin();
             index!= distributed_vector_short.End();
             ++index)
        {
            distributed_vector_short[index] = -0.5;
        }
        distributed_vector_short.Restore();

        DistributedVectorFactory factory2(2*number_nodes);
        Vec petsc_data_long = factory2.CreateVec();
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        for (DistributedVector::Iterator index = distributed_vector_long.Begin();
             index!= distributed_vector_long.End();
             ++index)
        {
            vm_stripe[index] = index.Global;
        }
        distributed_vector_long.Restore();

        writer.PutVector(ina_id, petsc_data_short);
        //Try to write striped data in the wrong columns
        TS_ASSERT_THROWS_ANYTHING(writer.PutStripedVector(vm_id, phi_e_id, petsc_data_long));
        //Try to write data of wrong size
        TS_ASSERT_THROWS_ANYTHING(writer.PutVector(ina_id, petsc_data_long));
        TS_ASSERT_THROWS_ANYTHING(writer.PutStripedVector(vm_id, ina_id, petsc_data_short));

        writer.Close();

        VecDestroy(petsc_data_long);
        VecDestroy(petsc_data_short);
    }

    /**
     * Tests copied (with some minor modifications) from TestColumnDataReaderWriter: to be refactored at some point
     */
    void TestDefineThings()
    {
        DistributedVectorFactory vec_factory(100);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time", "m secs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("T,i,m,e", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("", "msecs"));

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(6);

        // Data not increasing
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension(node_numbers, 100));
        node_numbers[2]=100;
        // Data is increasing but the last number is too large
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension(node_numbers, 100));

        mpTestWriter->DefineFixedDimension(5000);
        // Can't set fixed dimension more than once
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension(5000));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension(node_numbers, 100));

        int ina_var_id = 0;
        int ik_var_id = 0;
        int ik2_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineVariable("Dummy",""));

        // Defined twice
        TS_ASSERT_THROWS_ANYTHING(ik2_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));

        // Bad variable names/units
        TS_ASSERT_THROWS_ANYTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milli amperes"));
        TS_ASSERT_THROWS_ANYTHING(ik_var_id = mpTestWriter->DefineVariable("I   K", "milliamperes"));
        TS_ASSERT_THROWS_ANYTHING(ik_var_id = mpTestWriter->DefineVariable("I.K", "milliamperes"));
        TS_ASSERT_THROWS_ANYTHING(ik_var_id = mpTestWriter->DefineVariable("", "milliamperes"));

        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id, 1);

        delete mpTestWriter;
    }

    void TestEndDefineMode()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "testdefine"));

        // Ending define mode without having defined at least a variable and a fixed dimension should raise an exception
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->EndDefineMode());

        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));

        //In Hdf5 a fixed dimension should be defined always
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->EndDefineMode());
        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(60);
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension(node_numbers, 100));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineVariable("I_Ca", "milli amperes"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension(5000));

        // Can't call define fixed dimension again
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension(node_numbers, 100));

        // Test that we can't write incomplete data from a vector that doesn't have the right entries (0 to 59)
        DistributedVectorFactory factory(60);
        Vec petsc_data_short=factory.CreateVec();
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVector(0, petsc_data_short));
        VecDestroy(petsc_data_short);

        mpTestWriter->Close();
        delete mpTestWriter;
    }

    void TestCantAddUnlimitedAfterEndDefine()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "testdefine"));
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension(0));
        mpTestWriter->DefineFixedDimension(5000);

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        mpTestWriter->Close();
        delete mpTestWriter;
    }

    void TestAdvanceAlongUnlimitedDimension()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "testdefine"));

        int ina_var_id;
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension(5000));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutUnlimitedVariable(0.0));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->AdvanceAlongUnlimitedDimension());

        mpTestWriter->Close();
        delete mpTestWriter;
    }

    void TestCantWriteDataWhileInDefineMode()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory, "", "testdefine", false);
//        writer.DefineFixedDimension(number_nodes);
//
//        int node_id = writer.DefineVariable("Node","dimensionless");
//        int vm_id = writer.DefineVariable("V_m","millivolts");
//        int phi_e_id = writer.DefineVariable("Phi_e","millivolts");
//        int ina_id = writer.DefineVariable("I_Na","milliamperes");
//
//        writer.DefineUnlimitedDimension("Time", "msec");

//        writer.EndDefineMode();

        int node_id = 1;
        int vm_id = 2;
        int phi_e_id = 3;
        int ina_id = 4;

        DistributedVectorFactory factory(number_nodes);

        Vec petsc_data_short = vec_factory.CreateVec();
        DistributedVector distributed_vector_short = vec_factory.CreateDistributedVector(petsc_data_short);

        Vec node_number = vec_factory.CreateVec();
        DistributedVector distributed_node_number = vec_factory.CreateDistributedVector(node_number);

        for (DistributedVector::Iterator index = distributed_node_number.Begin();
             index!= distributed_node_number.End();
             ++index)
        {
            distributed_node_number[index] = index.Global;
            distributed_vector_short[index] = -0.5;
        }
        distributed_node_number.Restore();
        distributed_vector_short.Restore();

        DistributedVectorFactory factory2(2*number_nodes);
        Vec petsc_data_long = factory2.CreateVec();
        DistributedVector distributed_vector_long = factory2.CreateDistributedVector(petsc_data_long);

        for (DistributedVector::Iterator index = distributed_vector_long.Begin();
             index!= distributed_vector_long.End();
             ++index)
        {
            distributed_vector_long[index] =  1000 + index.Global;
        }
        distributed_vector_long.Restore();

        TS_ASSERT_THROWS_ANYTHING(writer.PutVector(node_id, node_number));
        TS_ASSERT_THROWS_ANYTHING(writer.PutVector(ina_id, petsc_data_short));
        TS_ASSERT_THROWS_ANYTHING(writer.PutStripedVector(vm_id, phi_e_id, petsc_data_long));
        TS_ASSERT_THROWS_ANYTHING(writer.PutUnlimitedVariable(0.0));
        TS_ASSERT_THROWS_ANYTHING(writer.AdvanceAlongUnlimitedDimension());

        writer.Close();
        VecDestroy(petsc_data_short);
        VecDestroy(node_number);
        VecDestroy(petsc_data_long);
    }

};
#endif /*TESTHDF5DATAWRITER_HPP_*/
