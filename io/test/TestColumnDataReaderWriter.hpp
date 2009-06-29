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


#ifndef _TESTCOLUMNDATAREADERWRITER_HPP_
#define _TESTCOLUMNDATAREADERWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "Exception.hpp"

using namespace std;

class TestColumnDataReaderWriter : public CxxTest::TestSuite
{

private:

    ColumnDataWriter *mpTestWriter;
    ColumnDataReader *mpTestReader;

    bool filesMatch(std::string testfileName, std::string goodfileName)
    {
        bool matching = true;

        ifstream testfile(testfileName.c_str(),ios::in);
        ifstream goodfile(goodfileName.c_str(),ios::in);
        std::string teststring;
        std::string goodstring;

        if (!testfile.is_open() || !goodfile.is_open())
        {
            TS_FAIL("Files not present.");
        }

        while (getline(testfile, teststring))
        {
            getline(goodfile,goodstring);
            if (teststring != goodstring)
            {
                matching = false;
            }
        }

        if (getline(goodfile,goodstring))
        {
            matching = false;
        }

        testfile.close();
        goodfile.close();
        return matching;
    }

public:

    void TestCreateColumnWriter()
    {
        // Create a new csvdata writer
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "test"));
        //check that the output directory exists
        //use the Boost libraries for this check
        //or stat()

        delete mpTestWriter;
    }

    void TestCreateColumnReader()
    {
        // File does not exist
        TS_ASSERT_THROWS_ANYTHING(mpTestReader = new ColumnDataReader("", "testdoesnotexist"));

        // File contains corrupt data
        TS_ASSERT_THROWS_ANYTHING(mpTestReader = new ColumnDataReader("io/test/data", "testbad", false));

        // .info file exists (unlimited) but _unlimited.dat file does not
        TS_ASSERT_THROWS_ANYTHING(mpTestReader = new ColumnDataReader("io/test/data", "UnlimitedMissing", false));

        // .info file exists (fixed dim) but .dat file does not
        TS_ASSERT_THROWS_ANYTHING(mpTestReader = new ColumnDataReader("io/test/data", "DatMissing", false));

        delete mpTestReader;
    }


    void TestDefineUnlimitedDimension()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time", "m secs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("T,i,m,e", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("", "msecs"));

        delete mpTestWriter;
    }

    void TestDefineFixedDimension()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless", 5000));

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension("Node ","dimensionless", 5000));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension("Node", "dimension.less", 5000));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension("*Node*","dimensionless", 5000));

        delete mpTestWriter;
    }

    void TestDefineVariable()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "test"));
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineVariable("Dummy", ""));

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
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "testdefine"));

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVariable(0, 0, 0));

        // Ending define mode without having defined a dimension and a variable should raise an exception
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->EndDefineMode());

        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineVariable("I_Ca", "milli amperes"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension("Node", "dimensionless", 5000));

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT(filesMatch(output_dir + "testdefine.dat", "io/test/data/testdefine_good.dat"));

        TS_ASSERT(filesMatch(output_dir + "testdefine.info", "io/test/data/testdefine_good.info"));
    }

    void TestCantAddUnlimitedAfterEndDefine()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "testdefine"));
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless", 0));
        mpTestWriter->DefineFixedDimension("Node","dimensionless", 5000);

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        delete mpTestWriter;
    }

    void TestPutVariableInUnlimitedFile()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "testunlimited"));
        int time_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        for (int i=0; i<=10; i++)
        {
            mpTestWriter->PutVariable(time_var_id, (double)(i)/10 - 0.1);
            mpTestWriter->PutVariable(ina_var_id, 12.0);
            mpTestWriter->PutVariable(ica_var_id, ((double)((i+1)*(i+1)))/3.0);
            mpTestWriter->PutVariable(ik_var_id, 7124.12355553*((double)(i+1))/12.0);
            mpTestWriter->AdvanceAlongUnlimitedDimension();
        }

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT_THROWS_NOTHING(mpTestReader = new ColumnDataReader("","testunlimited"));

        std::vector<double> values_ik = mpTestReader->GetValues("I_K");

        for (int i=0; i<11; i++)
        {
            TS_ASSERT_DELTA(values_ik[i]/(7124.12355553*((double)(i+1))/12.0), 1.0, 1e-3);
        }

        std::vector<double> time_values = mpTestReader->GetUnlimitedDimensionValues();
        for (int i=0; i < 10; i++)
        {
            TS_ASSERT_DELTA(time_values[i],i*0.1-0.1,1e-3);
        }

        // Test for coverage
        TS_ASSERT_THROWS_ANYTHING(values_ik = mpTestReader->GetValues("I_K", 3));
        TS_ASSERT_THROWS_ANYTHING(mpTestReader->GetValues("BadVar"));

        delete mpTestReader;
    }

    void TestPutNegativeVariable()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "testunlimitednegative"));
        int time_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_ANYTHING(time_var_id = mpTestWriter->DefineVariable("Time", "msecs"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.2));

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) -i));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 33.124));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 7124.12355553));
        mpTestWriter->AdvanceAlongUnlimitedDimension();
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -33.124));

        // Check that an incorrect var id causes an exception:
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVariable(234, -33.124));

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT(filesMatch(output_dir + "testunlimitednegative.dat", "io/test/data/testunlimitednegative_good.dat"));
    }

    void TestPutVariableInFixedFile()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "testfixed"));
        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        int short_id = 0;
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_ANYTHING(node_var_id = mpTestWriter->DefineVariable("Node", "dimensionless"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(short_id = mpTestWriter->DefineVariable("Short_column", "dimensionless"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVariable(node_var_id, 0, -1));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVariable(node_var_id, 0, -2));

        for (unsigned i=0; i<4; i++)
        {
            mpTestWriter->PutVariable(node_var_id, (double)(i+1), i);
            mpTestWriter->PutVariable(ina_var_id, 12.0, i);
            mpTestWriter->PutVariable(ica_var_id, ((double)((i+1)*(i+1)))/3.0, i);
            mpTestWriter->PutVariable(ik_var_id, 7124.12355553*((double)(i+1))/12.0, i);
        }

        for (unsigned i=0; i<2; i++)
        {
            mpTestWriter->PutVariable(short_id, (double)(i), i);
        }

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT(filesMatch(output_dir + "testfixed.dat", "io/test/data/testfixed_good.dat"));

        TS_ASSERT(filesMatch(output_dir + "testfixed.info", "io/test/data/testfixed_good.info"));

        TS_ASSERT_THROWS_NOTHING(mpTestReader = new ColumnDataReader("", "testfixed"));

        TS_ASSERT_THROWS_ANYTHING(mpTestReader->GetValues("BadVar", 0));

        for (int i=0; i<4; i++)
        {
            std::vector<double> values_ik = mpTestReader->GetValues("I_K",  i);
            TS_ASSERT_DELTA(values_ik[0]/(7124.12355553*((double)(i+1))/12.0), 1.0, 1e-3);
        }

        for (int i=0; i<4; i++)
        {
            std::vector<double> values_short = mpTestReader->GetValues("Short_column", i);
            if (i<2)
            {
                TS_ASSERT_DELTA(values_short[0], (double) i, 1e-3);
            }
            else
            {
                //Missing data in short column
                TS_ASSERT_DELTA(values_short[0], DBL_MAX, 1e-3);
            }
        }

        TS_ASSERT_THROWS_ANYTHING(std::vector<double> values_dodgy = mpTestReader->GetValues("non-existent_variable",1));

        // Check that get unlimited dimension values throws
        TS_ASSERT_THROWS_ANYTHING(std::vector<double> unlimited_values = mpTestReader->GetUnlimitedDimensionValues());

        delete mpTestReader;
    }

    void TestPutNegativeVariableInFixedFile()
    {

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "testfixed_negatives"));
        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_ANYTHING(node_var_id = mpTestWriter->DefineVariable("Node", "dimensionless"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        mpTestWriter->PutVariable(ina_var_id, (double) i,0);
        mpTestWriter->PutVariable(ina_var_id, (double) -i,1);
        mpTestWriter->PutVariable(ica_var_id, -33.124,3);
        mpTestWriter->PutVariable(ik_var_id, 7124.12355553,3);
        mpTestWriter->AdvanceAlongUnlimitedDimension();
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVariable(ica_var_id, -63.124,2));

        // Note: the above call to PutVariable will, in effect, execute AdvanceAlongUnlimitedDimension and
        //       therefore throw an exception, hence we have to repeat that call to PutVariable below
        mpTestWriter->PutVariable(ica_var_id, -63.124,2);
        mpTestWriter->PutVariable(node_var_id, 1,0);
        mpTestWriter->PutVariable(node_var_id, -4,3);

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT(filesMatch(output_dir + "testfixed_negatives.dat", "io/test/data/testfixed_negatives_good.dat"));
    }

    void TestPutVariableInFixedandUnlimitedFile()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("","testfixedandunlimited"));
        int time_var_id = 0;
        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_ANYTHING(time_var_id = mpTestWriter->DefineVariable("Time", "msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(node_var_id, 0, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i,0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i,1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -33.124,3));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 7124.12355553,3));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->AdvanceAlongUnlimitedDimension());

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(node_var_id, 0, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 63.124,2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -35.124,3));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.2));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVariable(time_var_id, 0.2,3));

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT(filesMatch(output_dir + "testfixedandunlimited_unlimited.dat",
                             "io/test/data/testfixedandunlimited_unlimited_good.dat"));

        TS_ASSERT(filesMatch(output_dir + "testfixedandunlimited_000001.dat",
                             "io/test/data/testfixedandunlimited_000001_good.dat"));

        TS_ASSERT_THROWS_NOTHING(mpTestReader = new ColumnDataReader("","testfixedandunlimited"));

        std::vector<double> time_values = mpTestReader->GetUnlimitedDimensionValues();
        std::vector<double> ica_values = mpTestReader->GetValues("I_Ca", 3);
        for (int i=0; i<2; i++)
        {
            TS_ASSERT_DELTA(time_values[i],(i+1)*0.1,1e-3);
            TS_ASSERT_DELTA(ica_values[i],-33.124 - i * 2,1e-3);
        }

        // Check exception thrown if dimension is not given
        TS_ASSERT_THROWS_ANYTHING(ica_values = mpTestReader->GetValues("I_Ca"));

        delete mpTestReader;
    }

    /*
     * This test is just to cover the line in ColumnDataWriter::PutVariable where
     * the fixed and unlimited dimensions are both set and the unlimited parameter
     * (i.e. time) is passed in negative.
     */
    void TestNegativeWithFixedAndUnlimitedDefined()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("", "testunlimitednegative2"));

        int time_var_id = 0;
        int node_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING( ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, -0.2));

        // Remember to delete - this closes the writer cleanly and means any data left
        // unwritten will be written to the datafile
        delete mpTestWriter;
    }

};

#endif //_TESTCOLUMNDATAREADERWRITER_HPP_
