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


#ifndef TEST1DCARDIACMECHANICSASSEMBLERS_HPP_
#define TEST1DCARDIACMECHANICSASSEMBLERS_HPP_

#include <cxxtest/TestSuite.h>
#include "Explicit1dCardiacMechanicsAssembler.hpp"
#include "Implicit1dCardiacMechanicsAssembler.hpp" 

class Test1dCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestUncoupledExplicit() throw(Exception)
    {
        Triangulation<1> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(7);
        
        Explicit1dCardiacMechanicsAssembler mechanics(&mesh);

        std::vector<double> active_tension(mechanics.GetTotalNumQuadPoints(), 0.5);
        
        mechanics.SetForcingQuantity( active_tension );
        mechanics.Solve(0, 0.1, 0.1); // times not used
        
        std::vector<Vector<double> > undeformed_position = mechanics.rGetUndeformedPosition();
        std::vector<Vector<double> > deformed_position = mechanics.rGetDeformedPosition();
        
        // not i=0 as X=0
        double factor = deformed_position[0](1)/undeformed_position[0](1);
        
        TS_ASSERT_LESS_THAN(factor, 1.0);
        
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            if(undeformed_position[0](i)!=0)
            {
                TS_ASSERT_DELTA(deformed_position[0](i)/undeformed_position[0](i), factor, 1e-4);
            }
        }
    }
    

    void TestUncoupledImplicit() throw(Exception)
    {
        Triangulation<1> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(7);
        
        Implicit1dCardiacMechanicsAssembler mechanics(&mesh);
        std::vector<double> caI(mechanics.GetTotalNumQuadPoints(), 0.02);
        
        mechanics.SetForcingQuantity( caI );
        mechanics.Solve(0, 0.1, 0.01);
        
        std::vector<Vector<double> > undeformed_position = mechanics.rGetUndeformedPosition();
        std::vector<Vector<double> > deformed_position = mechanics.rGetDeformedPosition();
        
        // not i=0 as X=0
        double factor = deformed_position[0](1)/undeformed_position[0](1);
        
        TS_ASSERT_LESS_THAN(factor, 1.0);
        
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            if(undeformed_position[0](i)!=0)
            {
                TS_ASSERT_DELTA(deformed_position[0](i)/undeformed_position[0](i), factor, 1e-4);
            }
        }
    }
    
};
#endif /*TEST1DCARDIACMECHANICSASSEMBLERS_HPP_*/
