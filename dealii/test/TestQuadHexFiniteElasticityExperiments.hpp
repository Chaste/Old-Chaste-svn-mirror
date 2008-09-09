/*

Copyright (C) University of Oxford, 2008

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

#ifndef TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_
#define TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_


#include <cxxtest/TestSuite.h>
#include "FiniteElasticityAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "FiniteElasticityTools.hpp"
#include "UblasCustomFunctions.hpp"


class TestQuadHexFiniteElasticityExperiments : public CxxTest::TestSuite
{
public:

    void TestConvergence() throw(Exception)
    {
        unsigned num_sims = 6;
        
        LogFile::Instance()->Set(1, "quad_hex_convergence");
        LogFile::Instance()->SetPrecision(9);
        
        unsigned num_elem_in_each_dir = 1;
        std::vector<unsigned> num_elems;
        std::vector<c_vector<double,2> > results;
        
        for(unsigned i=0; i<num_sims; i++)
        {
            Vector<double> body_force(2);
            body_force(0) = 0.06;
            MooneyRivlinMaterialLaw<2> mooney_rivlin_law(0.02);
    
            Triangulation<2> mesh;
            GridGenerator::hyper_cube(mesh, 0.0, 1.0);

            if(i>0)
            {
                mesh.refine_global(i);
            }

            FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
            
            std::stringstream dir;
            dir << "quad_hex_convergence/" << num_elem_in_each_dir;
    
            FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                           &mooney_rivlin_law,
                                                           body_force,
                                                           1.0,
                                                           dir.str());
            
            double tol = 1e-4*finite_elasticity.CalculateResidualNorm();
            if(tol > 1e-8)
            {
                tol = 1e-8;
            }
            else if(tol < 1e-12)
            {
                tol = 1e-12;
            }

            finite_elasticity.StaticSolve(true, tol);

            std::vector<Vector<double> >& r_solution = finite_elasticity.rGetDeformedPosition();
            std::vector<Vector<double> >& r_undef_position = finite_elasticity.rGetUndeformedPosition();

            unsigned top_corner_index = 2;
            TS_ASSERT_DELTA(r_undef_position[0](top_corner_index), 1.0, 1e-6);
            TS_ASSERT_DELTA(r_undef_position[1](top_corner_index), 1.0, 1e-6);

            LOG_AND_COUT(1,num_elem_in_each_dir << " " << r_solution[0](top_corner_index)<< " " << r_solution[1](top_corner_index));
            num_elems.push_back(num_elem_in_each_dir);
           
            c_vector<double,2> result;
            result(0) = r_solution[0](top_corner_index);
            result(1) = r_solution[1](top_corner_index);
            results.push_back(result);
            
            num_elem_in_each_dir*=2;
        }

        std::cout << "\n\n";
        for(unsigned i=0; i<results.size(); i++)
        {
            std::cout << std::setprecision(9) << num_elems[i] << " " << results[i](0) << " " << results[i](1) << std::endl;
        }
    }
};
#endif /*TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_*/
