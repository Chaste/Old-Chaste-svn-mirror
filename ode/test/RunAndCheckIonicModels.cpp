
#include <cxxtest/TestSuite.h>
#include <cmath>
#include "RunAndCheckIonicModels.hpp"

void RunOdeSolverWithIonicModel(AbstractCardiacCell *pOdeSystem,
                                double endTime,
                                std::string filename,
                                int stepPerRow,
                                bool doComputeExceptVoltage)
{
    double start_time = 0.0;
    
    if (doComputeExceptVoltage)
    {
        // Store the current system state
        std::vector<double> state_variables_ref = pOdeSystem->rGetStateVariables();
        std::vector<double> state_variables_copy = state_variables_ref;
        
        // Test ComputeExceptVoltage
        double v_init = pOdeSystem->GetVoltage();
        pOdeSystem->ComputeExceptVoltage(start_time, endTime);
        double v_end = pOdeSystem->GetVoltage();
        TS_ASSERT_DELTA(v_init, v_end, 1e-6);
        
        // Test SetVoltage
        pOdeSystem->SetVoltage(1e6);
        TS_ASSERT_DELTA(pOdeSystem->GetVoltage(), 1e6, 1e-6);
        
        // Reset the system
        pOdeSystem->SetStateVariables(state_variables_copy);
    }
    
    // Solve
    OdeSolution solution = pOdeSystem->Compute(start_time, endTime);
    
    // Write data to a file using ColumnDataWriter
    SaveSolution(filename, pOdeSystem, solution, stepPerRow);
}

void SaveSolution(std::string baseResultsFilename, AbstractCardiacCell *pOdeSystem,
                  OdeSolution& rSolution, int stepPerRow)
{
    // Write data to a file using ColumnDataWriter
    
    ColumnDataWriter writer("TestIonicModels",baseResultsFilename,false);
    int time_var_id = writer.DefineUnlimitedDimension("Time","ms");
    
    std::vector<int> var_ids;
    for (unsigned i=0; i<pOdeSystem->rGetVariableNames().size(); i++)
    {
        var_ids.push_back(writer.DefineVariable(pOdeSystem->rGetVariableNames()[i],
                                                pOdeSystem->rGetVariableUnits()[i]));
    }
    writer.EndDefineMode();
    
    for (unsigned i = 0; i < rSolution.rGetSolutions().size(); i+=stepPerRow)
    {
        writer.PutVariable(time_var_id, rSolution.rGetTimes()[i]);
        for (unsigned j=0; j<var_ids.size(); j++)
        {
            writer.PutVariable(var_ids[j], rSolution.rGetSolutions()[i][j]);
        }
        writer.AdvanceAlongUnlimitedDimension();
    }
    writer.Close();
}

void CheckCellModelResults(std::string baseResultsFilename)
{
    /*
     * Check the cell model against a previous version
     * or another source e.g. Alan's COR
     */
    
    // read data entries for the new file and compare to valid data from
    // other source
    ColumnDataReader data_reader("TestIonicModels", baseResultsFilename);
    std::vector<double> times = data_reader.GetValues("Time");
    std::vector<double> voltages = data_reader.GetValues("V");
    ColumnDataReader valid_reader("ode/test/data", baseResultsFilename+"ValidData",
                                  false);
    std::vector<double> valid_times = valid_reader.GetValues("Time");
    std::vector<double> valid_voltages = valid_reader.GetValues("V");
    
    TS_ASSERT_EQUALS(times.size(), valid_times.size());
    for (unsigned i=0; i<valid_times.size(); i++)
    {
        TS_ASSERT_DELTA(times[i], valid_times[i], 1e-12);
        // adjust tol to data
        TS_ASSERT_DELTA(voltages[i], valid_voltages[i], 1e-6);
    }
}

void CompareCellModelResults(std::string baseResultsFilename1, std::string baseResultsFilename2, double tolerance)
{
    // Compare 2 sets of results, e.g. from 2 different solvers for the same model.
    // Initially we assume the time series are the same; this will change.
    // If the time series differ, the finer resolution must be given first.
    ColumnDataReader data_reader1("TestIonicModels", baseResultsFilename1);
    std::vector<double> times1 = data_reader1.GetValues("Time");
    std::vector<double> voltages1 = data_reader1.GetValues("V");
    std::vector<double> calcium1 = data_reader1.GetValues("CaI");
    std::vector<double> h1 = data_reader1.GetValues("h");
    
    ColumnDataReader data_reader2("TestIonicModels", baseResultsFilename2);
    std::vector<double> times2 = data_reader2.GetValues("Time");
    std::vector<double> voltages2 = data_reader2.GetValues("V");
    std::vector<double> calcium2 = data_reader2.GetValues("CaI");
    std::vector<double> h2 = data_reader2.GetValues("h");
    
    TS_ASSERT(times1.size() >= times2.size());
    double last_v = voltages2[0];
    double tol = tolerance;
    for (unsigned i=0, j=0; i<times2.size(); i++)
    {
        // Find corresponding time index
        while (j<times1.size() && times1[j] < times2[i] - 1e-12)
        {
            j++;
        }
        
        // Set tolerance higher in upstroke
        if (fabs(voltages2[i] - last_v) > 0.05)
        {
            tol = tolerance * 25;
        }
        else
        {
            tol = tolerance;
        }
        last_v = voltages2[i];
        
        TS_ASSERT_DELTA(times1[j], times2[i], 1e-12);
        // adjust tol to data
        TS_ASSERT_DELTA(voltages1[j], voltages2[i], tol);
        TS_ASSERT_DELTA(calcium1[j],  calcium2[i],  tol/1000);
        TS_ASSERT_DELTA(h1[j],        h2[i],        tol/10);
    }
}

