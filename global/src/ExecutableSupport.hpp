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

#ifndef EXECUTABLESUPPORT_HPP_
#define EXECUTABLESUPPORT_HPP_

/**
 * Various helpful static methods for people writing their own executables
 * within the Chaste framework.
 *
 * Most executables will just need to call StandardStartup as the first thing in their
 * main() function, and FinalizePetsc before quitting.  The other methods allow you to
 * fine-tune what output is presented to users.
 */
class ExecutableSupport
{
public:
    /**
     * Initialise PETSc from the command line arguments.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void InitializePetsc(int* pArgc, char*** pArgv);

    /**
     * Display Chaste's copyright information.
     */
    static void ShowCopyright();

    /**
     * Output extra diagnostics when Chaste is launched in parallel.
     */
    static void ShowParallelLaunching();

    /**
     * Call InitializePetsc, ShowCopyright, then ShowParallelLaunching.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void StandardStartup(int* pArgc, char*** pArgv);

    /**
     * Shut down PETSc so we exit cleanly.
     */
    static void FinalizePetsc();
};

#endif /* EXECUTABLESUPPORT_HPP_ */
