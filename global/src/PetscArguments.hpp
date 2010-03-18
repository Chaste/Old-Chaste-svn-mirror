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

#ifndef PETSCARGUMENTS_HPP_
#define PETSCARGUMENTS_HPP_

/**
 * A convenient holder for the command line arguments.
 *
 * The cxxtest harness will fill in the member variables when a test is
 * started.  They can then be read by PETSc when it is initialised.
 */
class PetscArguments
{
public:

    /** The number of command line arguments. */
    int* p_argc;

    /** The arguments themselves. */
    char*** p_argv;

    /** Get the single instance of this class. */
    static PetscArguments* Instance();

private:

    /** Default constructor. Should never be called directly, call PetscArguments::Instance() instead.*/
    PetscArguments();

    /** Copy constructor. */
    PetscArguments(const PetscArguments&);

    /** Overloaded assignment operator. */
    PetscArguments& operator= (const PetscArguments&);

    /** The single instance of the class. */
    static PetscArguments* mpInstance;
};

#endif // PETSCARGUMENTS_HPP_
