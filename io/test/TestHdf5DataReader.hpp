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
#ifndef TESTHDF5READER_HPP_
#define TESTHDF5READER_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "Hdf5DataWriter.hpp"
#include "Hdf5DataReader.hpp"
#include "DistributedVectorFactory.hpp"

class TestHdf5DataReader : public CxxTest::TestSuite
{
private :

    std::string h5file_name;
    #define DATASETNAME "IntArray"

    void WriteDataTestSimpleReadDirectlyWithHdf5()
    {
        int const NX = 5;                      /* dataset dimensions */
        int const NY = 6;
        int const RANK = 2;

        h5file_name = OutputFileHandler::GetChasteTestOutputDirectory() + "SDS.h5";

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
        file = H5Fcreate(h5file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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
        dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT);

        /*
         * Write the data to the dataset using default transfer properties.
         */
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

        /*
         * Close/release resources.
         */
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Fclose(file);
        PetscTools::Barrier("WriteDataTestSimpleReadDirectlyWithHdf5");
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
                for (k = 0; k < NZ; k++)
                    data_out[j][i][k] = 0;
            }
        }

        // Open the file and the dataset
        file = H5Fopen(h5file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
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

        // Define hyperslab in the dataset
        offset[0] = 1;
        offset[1] = 2;
        count[0]  = NX_SUB;
        count[1]  = NY_SUB;
        status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
                     count, NULL);

        // Define the memory dataspace
        dimsm[0] = NX;
        dimsm[1] = NY;
        dimsm[2] = NZ;
        memspace = H5Screate_simple(RANK_OUT,dimsm,NULL);

        // Define memory hyperslab
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

        // Close/release resources
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Sclose(memspace);
        H5Fclose(file);
    }

private:

    static const unsigned NUMBER_NODES = 100;

    void WriteMultiStepData()
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataWriter writer(factory, "hdf5_reader", "hdf5_test_complete_format", false);
        writer.DefineFixedDimension(NUMBER_NODES);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int ik_id = writer.DefineVariable("I_K", "milliamperes");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        Vec petsc_data_1=factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2=factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3=factory.CreateVec();
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

            // Don't advance in the last iteration
            if (time_step < 9)
            {
                writer.AdvanceAlongUnlimitedDimension();
            }
        }
        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);
        VecDestroy(petsc_data_3);
        writer.Close();
    }

public:

    void TestMultiStepReader() throw (Exception)
    {
        WriteMultiStepData();

        Hdf5DataReader reader("hdf5_reader", "hdf5_test_complete_format");

        std::vector<std::string> variable_names=reader.GetVariableNames();
        TS_ASSERT_EQUALS(variable_names.size(), 3U);
        TS_ASSERT_EQUALS(variable_names[0], "Node");
        TS_ASSERT_EQUALS(reader.GetUnit("Node"), "dimensionless");
        TS_ASSERT_EQUALS(variable_names[1], "I_K");
        TS_ASSERT_EQUALS(reader.GetUnit("I_K"), "milliamperes");
        TS_ASSERT_EQUALS(variable_names[2], "I_Na");
        TS_ASSERT_EQUALS(reader.GetUnit("I_Na"), "milliamperes");

        for (unsigned node_index=0; node_index<NUMBER_NODES; node_index++)
        {
            std::vector<double> node_values = reader.GetVariableOverTime("Node", node_index);
            std::vector<double> i_k_values = reader.GetVariableOverTime("I_K", node_index);
            std::vector<double> i_na_values = reader.GetVariableOverTime("I_Na", node_index);

            TS_ASSERT_EQUALS(node_values.size(), 10u);
            TS_ASSERT_EQUALS(i_k_values.size(), 10u);
            TS_ASSERT_EQUALS(i_na_values.size(), 10u);

            for (unsigned i=0; i<node_values.size(); i++)
            {
                TS_ASSERT_DELTA( node_values[i], node_index, 1e-9);
                TS_ASSERT_DELTA( i_k_values[i], i*1000 + 100 + node_index, 1e-9);
                TS_ASSERT_DELTA( i_na_values[i], i*1000 + 200 + node_index, 1e-9);
            }
        }

        unsigned NUMBER_NODES=100;
        DistributedVectorFactory factory(NUMBER_NODES);

        Vec petsc_data_1=factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2=factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3=factory.CreateVec();
        DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_3);

        TS_ASSERT_EQUALS(reader.GetNumberOfRows(), NUMBER_NODES);
        for (unsigned time_step=0; time_step<10; time_step++)
        {
            reader.GetVariableOverNodes(petsc_data_1, "Node", time_step);
            reader.GetVariableOverNodes(petsc_data_2, "I_K", time_step);
            reader.GetVariableOverNodes(petsc_data_3, "I_Na", time_step);

            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();

            // Check values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index!= distributed_vector_1.End();
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

        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);
        VecDestroy(petsc_data_3);
        reader.Close();
    }

    void TestNonMultiStepExceptions ()
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataWriter writer(factory, "hdf5_reader", "hdf5_test_overtime_exceptions", false);
        writer.DefineFixedDimension(NUMBER_NODES);

        writer.DefineVariable("Node", "dimensionless");
        writer.DefineVariable("I_K", "milliamperes");
        writer.DefineVariable("I_Na", "milliamperes");

        writer.EndDefineMode();
        writer.Close();

        Hdf5DataReader reader("hdf5_reader", "hdf5_test_overtime_exceptions");

        // Unlimited dimension has a default return value
        std::vector<double> times = reader.GetUnlimitedDimensionValues();
        TS_ASSERT_EQUALS(times.size(),1U);
        TS_ASSERT_EQUALS(times[0],0.0);

        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("Node", 99/*node*/),
                "The file does not contain time dependent data");

        Vec data = factory.CreateVec();
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "Node", 1/*timestep*/),
                "The file does not contain time dependent data");

        VecDestroy(data);
        reader.Close();
        
        // Missing dir
        FileFinder absent_dir("absent_dir", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(Hdf5DataReader(absent_dir, "base"), "Directory does not exist: ");
    }

    void TestMultiStepExceptions () throw (Exception)
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataWriter writer(factory, "hdf5_reader", "hdf5_test_overtime_exceptions", false);
        DistributedVectorFactory vec_factor(NUMBER_NODES);
        writer.DefineFixedDimension(NUMBER_NODES);

        writer.DefineVariable("Node", "dimensionless");
        writer.DefineVariable("I_K", "milliamperes");
        writer.DefineVariable("I_Na", "milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        writer.AdvanceAlongUnlimitedDimension();

        writer.Close();
        TS_ASSERT_THROWS_CONTAINS(Hdf5DataReader reader2("hdf5_reader", "hdf5_wrong_name"),"Hdf5DataReader could not open ");
        Hdf5DataReader reader("hdf5_reader", "hdf5_test_overtime_exceptions");

        TS_ASSERT_THROWS_NOTHING(reader.GetVariableOverTime("Node", 99/*node*/));
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("WrongName", 99/*node*/),
                "The file doesn\'t contain data for variable WrongName");
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("Node", 100/*node*/),
                "The file doesn\'t contain info of node 100");

        Vec data = factory.CreateVec();
        reader.GetVariableOverNodes(data, "Node", 0/*timestep*/);
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "WrongName"),
                "The file does not contain data for variable WrongName"); //Wrong name
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "I_K", 1/*timestep*/),
                "The file does not contain data for timestep number 1"); //Time step doesn't exist

        DistributedVectorFactory factory2(NUMBER_NODES+1);
        Vec data_too_big = factory2.CreateVec();
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data_too_big, "Node", 0/*timestep*/),
                "Could not read data because Vec is the wrong size"); //Data too big

        VecDestroy(data);
        VecDestroy(data_too_big);
        reader.Close();
    }


    void TestIncompleteData() throw (Exception)
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataReader reader("io/test/data","hdf5_test_full_format_incomplete", false);

        std::vector<std::string> variable_names = reader.GetVariableNames();
        TS_ASSERT_EQUALS(variable_names.size(), 3u);
        TS_ASSERT_EQUALS(variable_names[0], "Node");
        TS_ASSERT_EQUALS(reader.GetUnit("Node"), "dimensionless");
        TS_ASSERT_EQUALS(variable_names[1], "I_K");
        TS_ASSERT_EQUALS(reader.GetUnit("I_K"), "milliamperes");
        TS_ASSERT_EQUALS(variable_names[2], "I_Na");
        TS_ASSERT_EQUALS(reader.GetUnit("I_Na"), "milliamperes");

        // Can't read into a PETSc Vec
        Vec data = factory.CreateVec();
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "Node", 1/*timestep*/),
                "You can only get a vector for complete data");
        VecDestroy(data);

        std::vector<unsigned> nodes=reader.GetIncompleteNodeMap();
        TS_ASSERT_EQUALS(nodes.size(), 3U);
        TS_ASSERT_EQUALS(nodes[0], 21U);
        TS_ASSERT_EQUALS(nodes[1], 47U);
        TS_ASSERT_EQUALS(nodes[2], 60U);

        // Can read one of the nodes that was written
        std::vector<double> twenty_one = reader.GetVariableOverTime("Node", 21);
        TS_ASSERT_EQUALS(twenty_one.size(), 10u);
        for (unsigned i=0; i<twenty_one.size(); i++)
        {
            TS_ASSERT_EQUALS(twenty_one[i], 21u);
        }
        // Can read more of the data
        std::vector<double> na_47 = reader.GetVariableOverTime("I_Na", 47);
        TS_ASSERT_EQUALS(na_47.size(), 10u);
        for (unsigned i=0; i<na_47.size(); i++)
        {
            TS_ASSERT_EQUALS(na_47[i], i*1000u + 200u + 47u);
        }

        // Data not included
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("Node", 22),"The incomplete file does not contain info of node 22");
    }
};

#endif /*TESTHDF5READER_HPP_*/
