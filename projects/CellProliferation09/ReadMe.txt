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

CHASTE CELLPROLIFERATION09 PROJECT FOLDER

A collection of test script files and their results, which 
use the central chaste cancer code to recreate the results 
presented in our J. Cell Proilferation (2009) paper.

There are three folders - build, src and test. 

The <build> folder will contain the executables that you compile and can be ignored.

The <src> folder contains:

	* Hello.hpp, Hello.cpp - an example class, used by TestHello.hpp
	
    * SunterSetup.hpp, SunterSetup.cpp - helper class that sets up a crypt as described by Sunter et al. (1979) 
				        				 for different sites along the length of the colon

The <test> folder contains:

	* TestHello.hpp - run TestHello.hpp to check your Chaste-CellProliferation09 installation. 
					  The rest of the necessary source files are in the main chaste folder.

    * TestGeneratePlotsOfASingleCell.hpp - This file tracks the progeny of a single cell in the crypt over several experiments.

	* TestGenerateSteadyStateCryptCellProliferation.hpp - This file can be run with 
			different options to set up the steady-state crypt archives, which the main simulations are then run from.
			
 	* TestMeinekeLabellingExperimentsSunterData.hpp - This file simulates the crypt labelling and sections experiments.
 	
 	* TestMutationSpread.hpp - A single monoclonality experiment, this generates results to make movie files.
 	
 	* TestNicheSuccessionTimeDistributions.hpp - Runs many monoclonality experiments to generate succession time data.
 	
 	* <data> folder - contains the various steady-state crypt archives.
 	
 	* <results> folder - contains all of the data (produced by above test files), scripts and matlab files 
 						 necessary to generate the figures in the CellProliferation09 paper. 


To plot results using the .m files in the test/results MatLab(TM) you will first need to add the chaste matlab folder 
to your matlab path, using a MatLab command of the form:
addpath(<YOUR PATH>/chaste/anim/matlab);





