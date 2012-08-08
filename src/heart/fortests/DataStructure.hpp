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

#ifndef DATASTRUCTURE_HPP_
#define DATASTRUCTURE_HPP_

#include <fstream>
#include <vector>

#include "UblasVectorInclude.hpp"
#include "Exception.hpp"
#include "FileFinder.hpp"

/**
 * Helper class to read in data
 */
class DataStructure
{
  private:
    std::vector<std::string> mDrugNames;
    std::vector<unsigned> mRedfernCategory;
    std::vector<c_vector<double,3> > mIc50values;
    std::vector<c_vector<double,2> > mClinicalDoseRange;
    std::vector<double> mGrandiMeasure;

  public:
    DataStructure(std::string fileName)
    {
        LoadDrugDataFromFile(fileName);
    }

    DataStructure(FileFinder& rFileFinder)
    {
        LoadDrugDataFromFile(rFileFinder.GetAbsolutePath());
    }

    std::string GetDrugName(unsigned drugIndex)
    {
        assert(drugIndex < GetNumDrugs());
        return mDrugNames[drugIndex];
    }

    unsigned GetRedfernCategory(unsigned drugIndex)
    {
        assert(drugIndex < GetNumDrugs());
        return mRedfernCategory[drugIndex];
    }

    double GetGrandiMeasure(unsigned drugIndex)
    {
        assert(drugIndex < GetNumDrugs());
        if (mGrandiMeasure[drugIndex] < -998)
        {
            EXCEPTION("No data available on Grandi measure for this drug.");
        }
        return mGrandiMeasure[drugIndex];
    }

    double GetIC50Value(unsigned drugIndex, unsigned channelIndex)
    {
        assert(drugIndex < GetNumDrugs());
        assert(channelIndex<3u);
        double ic50 = mIc50values[drugIndex](channelIndex);

        /**
         * If an IC50 value takes the value
         *  -1 in the data file that means "affect unknown".
         *  -2 in the data file means "known to have no affect".
         */
        if (fabs(ic50 + 2) < 1e-9)
        {   // This drug is known to have no effect on this channel
            // We set the IC50 value to DBL_MAX for the block calculations.
            ic50 = DBL_MAX;
        }

        return ic50;
    }

    double GetClinicalDoseRange(unsigned drugIndex, unsigned lowOrHigh)
    {
        assert(lowOrHigh==0 || lowOrHigh==1);
        assert(drugIndex < GetNumDrugs());
        if (mClinicalDoseRange[drugIndex](0) < 0)
        {
            EXCEPTION("No data available on clinical dose for this drug.");
        }
        return mClinicalDoseRange[drugIndex](lowOrHigh);
    }

    void LoadDrugDataFromFile(std::string fileName)
    {

        std::ifstream indata; // indata is like cin
        indata.open(fileName.c_str()); // opens the file
        if(!indata)
        { // file couldn't be opened
            EXCEPTION("Couldn't open data file: " + fileName);
        }

        bool line_is_not_blank = true;
        unsigned line_counter = 0;
        while (line_is_not_blank)
        {
            std::string this_line;
            getline(indata, this_line);
            std::stringstream line(this_line);

            unsigned category;
            std::string name;
            std::string in_redfern_figs;
            c_vector<double, 3> ic50s;
            c_vector<double, 2> doses;
            double grandi_measure;
            ic50s.clear();
            doses.clear();

            line >> name;
            line >> category;
            line >> ic50s(0);
            line >> ic50s(1);
            line >> ic50s(2);
            line >> doses(0);
            line >> doses(1);
            line >> in_redfern_figs;
            line >> grandi_measure;

            mDrugNames.push_back(name);
            mRedfernCategory.push_back(category);
            mIc50values.push_back(ic50s);
            mClinicalDoseRange.push_back(doses);
            mGrandiMeasure.push_back(grandi_measure);

            line_counter++;
            //std::cout << "line_counter = " <<line_counter << "\n" << std::flush;
            if (indata.eof())
            {
                break;
            }
        }

    }

    unsigned GetNumDrugs(void)
    {
        assert(mDrugNames.size()==mRedfernCategory.size());
        assert(mDrugNames.size()==mIc50values.size());
        assert(mDrugNames.size()==mClinicalDoseRange.size());
        return mDrugNames.size();
    }

    /**
     * Calculate the probability of a channel being open given this drug, IC50 and hill coefficient.
     *
     * Note: A negative IC50 value is interpreted as "drug has no effect on this channel".
     *
     * @param rConc  concentration of the drug.
     * @param rIC50  IC50 value for this drug and channel
     * @param hill  Hill coefficient for this drug dependent inactivation curve (defaults to 1).
     *
     * @return proportion of channels which are still active at this drug concentration
     */
    static double CalculateConductanceFactor(const double& rConc, const double& rIC50, double hill = 1.0)
    {
        if (rIC50 < 0)
        {
            return 1.0;
        }
        else
        {
            return 1.0/(1.0 + pow((rConc/rIC50), hill));
        }
    }
};

#endif // DATASTRUCTURE_HPP_
