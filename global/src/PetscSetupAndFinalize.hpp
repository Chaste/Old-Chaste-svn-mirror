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


#ifndef _PETSCSETUPANDFINALIZE_HPP_
#define _PETSCSETUPANDFINALIZE_HPP_

/**
 * This file is designed to be included by any test suites that use PETSc.
 * It controls the PETSc initialisation and finalisation.
 *
 * Currently it will dump info on any non-freed vectors or matrices
 * on finalisation.
 */

#include <cxxtest/GlobalFixture.h>
#include <petsc.h>
#include <stdlib.h>

#include <unistd.h>
#include <iostream>

#include "PetscException.hpp"

class PetSCSetup : public CxxTest::GlobalFixture
{
public:
    /// Standard setup method for PETSc.
    bool setUpWorld()
    {
        /**
         * The cxxtest_argc and cxxtest_argv variables are global, and filled in
         * from the arguments passed to the cxxtest test suite runner.
         */
        PETSCEXCEPT(PetscInitialize(&cxxtest_argc, &cxxtest_argv,
                                    PETSC_NULL, PETSC_NULL) );

        char buf[10000];
        std::cout << std::endl << "CWD: " << getcwd(buf, 10000) << std::endl;
        std::cout << "Root: " << CHASTE_ROOT << std::endl;
        chdir(CHASTE_ROOT);
        std::cout << "CWD: " << getcwd(buf, 10000) << std::endl;

        return true;
    }
    // Clean up PETSc after running all tests.
    bool tearDownWorld()
    {
        //std::cout << "Finalizing..." << std::endl;
        PETSCEXCEPT(PetscFinalize());
        return true;
    }
};
static PetSCSetup thisSetup;


#endif //_PETSCSETUPANDFINALIZE_HPP_
