/*

Copyright (C) University of Oxford, 2005-2011

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
#ifndef TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_
#define TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_


#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ToyCompressibleMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "NonlinearElasticityTools.hpp"





/// todos: #1699:
///   more compressibility tests, nonlinear (I3) material law, (possible material law refactor)
///   stop using linear systems, two matrices in compressible case, matrix memory allocation




class TestCompressibleNonlinearElasticitySolver : public CxxTest::TestSuite
{
public:
    // This is purely for coverage of assembling a 3D system (and also uses alternative, heterogeneous
    // constructor, also for coverage)
    void TestAssembleSystem3D() throw (Exception)
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);

        mesh.ConstructFromMeshReader(mesh_reader1);

        ToyCompressibleMaterialLaw<3> law(1.0, 0.0, -1.0);
        std::vector<AbstractMaterialLaw<3>*> laws;
        laws.push_back(&law);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        CompressibleNonlinearElasticitySolver<3> solver(&mesh,
                                                        laws,
                                                        zero_vector<double>(3),
                                                        1.0,
                                                        "",
                                                        fixed_nodes);


        solver.AssembleSystem(true, true);
    }

    // It just tests that nothing happens if zero force and tractions are given
    void TestWithZeroDisplacement() throw(Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        ToyCompressibleMaterialLaw<2> law(1.0, 0.0, -1.0);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        CompressibleNonlinearElasticitySolver<2> solver(&mesh,
                                                        &law,
                                                        zero_vector<double>(2),
                                                        1.0,
                                                        "",
                                                        fixed_nodes);

        solver.Solve();
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u);

        // get deformed position
        std::vector<c_vector<double,2> >& r_deformed_position
            = solver.rGetDeformedPosition();

        for (unsigned i=0; i<r_deformed_position.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], r_deformed_position[i](0), 1e-8);
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[1], r_deformed_position[i](1), 1e-8);
        }
    }

    /**
     *  Suppose the deformation is given to be x = (alpha X, beta Y), and the material law is the toy
     *  law W(I1,I2,I3) = c1(I1-3) + c3(I3-1)   (where c3=-c1).
     *  On the unit square we specify displacement boundaries on the X=0 which match the deformation,
     *  then the deformation which corresponds to constant F=diag(a,b) => constant 1st PK stress
     *  S=diag(2*c1*alpha + 2*c3/alpha, 2*c1*beta + 2*c3/beta) if we choose
     *   - zero body force
     *   - beta = sqrt(-c3/c1)  (=> S22=0) and zero traction bcs on lower/upper surfaces
     *   - a traction bc of 2*c1*alpha + 2*c3/alpha (=S11) on the X=1 surface
     */
    void TestSolveForSimpleDeformation() throw(Exception)
    {
        double c1 = 1.0;
        double c3 = -1.0;
        double alpha = 0.9;
        double beta = sqrt(-c3/c1);
        double traction_value = 2*c1*alpha + 2*c3/alpha;

        c_vector<double,2> body_force = zero_vector<double>(2);
        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        ToyCompressibleMaterialLaw<2> law(c1, 0.0, c3);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
                c_vector<double,2> new_position;
                new_position(0) = 0;
                new_position(1) = beta*mesh.GetNode(i)->rGetLocation()[1];
                locations.push_back(new_position);
            }
        }

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = traction_value;
        traction(1) = 0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        assert(boundary_elems.size()==num_elem);

        CompressibleNonlinearElasticitySolver<2> solver(&mesh,
                                                        &law,
                                                        body_force,
                                                        1.0,
                                                        "comp_nonlin_elas_non_zero_bcs",
                                                        fixed_nodes,
                                                        &locations);

        solver.SetSurfaceTractionBoundaryConditions(boundary_elems, tractions);

        // coverage
        solver.SetKspAbsoluteTolerance(1e-10);

        solver.Solve();

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 3u); // 'hardcoded' answer, protects against jacobian getting messed up

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        for (unsigned i=0; i<fixed_nodes.size(); i++)
        {
            unsigned index = fixed_nodes[i];
            TS_ASSERT_DELTA(r_solution[index](0), locations[i](0), 1e-8);
            TS_ASSERT_DELTA(r_solution[index](1), locations[i](1), 1e-8);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = alpha*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = beta*mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
        }

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }
};

#endif /* TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_ */
