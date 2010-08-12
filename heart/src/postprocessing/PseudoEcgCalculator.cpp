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

#include "PseudoEcgCalculator.hpp"
#include "HeartConfig.hpp"
#include "PetscTools.hpp"
#include "Version.hpp"
#include <iostream>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> ::GetIntegrand(ChastePoint<SPACE_DIM>& rX,
                                c_vector<double,PROBLEM_DIM>& rU,
                                c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU)
{
    c_vector<double,SPACE_DIM> r_vector = rX.rGetLocation()- mProbeElectrode.rGetLocation();
    double norm_r = norm_2(r_vector);
    if (norm_r <= DBL_EPSILON)
    {
        EXCEPTION("Probe is on a mesh Gauss point.");
    }
    c_vector<double,SPACE_DIM> grad_one_over_r = - (r_vector)*SmallPow( (1/norm_r) , 3);
    matrix_row<c_matrix<double, PROBLEM_DIM, SPACE_DIM> > grad_u_row(rGradU, 0);
    double integrand = inner_prod(grad_u_row, grad_one_over_r);

    return -mDiffusionCoefficient*integrand;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> ::PseudoEcgCalculator (AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                                                                                 const ChastePoint<SPACE_DIM>& rProbeElectrode,
                                                                                 std::string directory,
                                                                                 std::string hdf5File,
                                                                                 std::string variableName,
                                                                                 bool makeAbsolute)
                                      : mrMesh(rMesh),
                                        mProbeElectrode(rProbeElectrode),
                                        mVariableName(variableName)

{
    
    mpDataReader = new Hdf5DataReader(directory, hdf5File, makeAbsolute);
    mNumberOfNodes = mpDataReader->GetNumberOfRows();
    mNumTimeSteps = mpDataReader->GetVariableOverTime(mVariableName, 0u).size();
    mDiffusionCoefficient = 1.0;
    //check that the hdf file was generated by simulations from the same mesh
    assert(mNumberOfNodes == mrMesh.GetNumNodes());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~PseudoEcgCalculator()
{
    delete mpDataReader;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetDiffusionCoefficient(double diffusionCoefficient)
{
    assert(diffusionCoefficient>=0);
    mDiffusionCoefficient = diffusionCoefficient;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputePseudoEcgAtOneTimeStep (unsigned timeStep)
{
    Vec solution_at_one_time_step = PetscTools::CreateVec(mNumberOfNodes);
    mpDataReader->GetVariableOverNodes(solution_at_one_time_step, mVariableName , timeStep);

    double pseudo_ecg_at_one_timestep;
    try
    {
        pseudo_ecg_at_one_timestep = Calculate(mrMesh, solution_at_one_time_step);
    }
    catch (Exception &e)
    {
        VecDestroy(solution_at_one_time_step);
        throw e;
    }
    VecDestroy(solution_at_one_time_step);
    return pseudo_ecg_at_one_timestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::WritePseudoEcg ()
{
    std::stringstream stream;
    stream << "PseudoEcg.dat";
    OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
    out_stream p_file=out_stream(NULL);
    if (PetscTools::AmMaster())
    {
        p_file = output_file_handler.OpenOutputFile(stream.str());
    }
    for (unsigned i = 0; i < mNumTimeSteps; i++)
    {
        double pseudo_ecg_at_one_timestep = ComputePseudoEcgAtOneTimeStep(i);
        if (PetscTools::AmMaster())
        {
            *p_file << pseudo_ecg_at_one_timestep << "\n";
        }
    }
    if (PetscTools::AmMaster())
    {
        //write provenance info
	    std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();
	    *p_file << comment;
        p_file->close();
    }
}
/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class PseudoEcgCalculator<1,1,1>;
//template class PseudoEcgCalculator<1,2,1>;
//template class PseudoEcgCalculator<1,3,1>;
template class PseudoEcgCalculator<2,2,1>;
//template class PseudoEcgCalculator<2,3,1>;
template class PseudoEcgCalculator<3,3,1>;

