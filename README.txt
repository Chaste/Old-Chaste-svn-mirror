This folder contains all the code used to generate the results presented in

Mirams et al. (2011), 
Simulation of multiple ion channel block provides improved early prediction of compounds' clinical torsadogenic risk.
Cardiovascular Research, doi:10.1093/cvr/CVR044

It is packaged as a "bolt-on" project for Chaste, compatible with Chaste version 2.1 only.

You should first download the Chaste v2.1 source code from http://www.comlab.ox.ac.uk/chaste/download.html, 
and then unpack this folder into Chaste/projects/CardiovascRes11

(It is important to use this foldername so that the drug data file can be found)

The code can then be compiled and tested from the Chaste source folder using

scons cl=1 projects/CardiovascRes11

To recreate the paper results then run:
build/<BUILD_TYPE>/TestModifyingConductancesRunner
with the appropriate options (it will prompt for these).

Multiple runs of TestModifyingConductancesRunnerruns then compile output into
$CHASTE_TEST_OUTPUT/VaryingConductances

Matlab scripts are then used to perform the postprocessing for e.g. restitution slopes,
and these will store those results in 
test/matlab_analysis/results

Classification is performed on those results by scripts in
test/matlab_analysis/classification

In addition there is an executable TorsadePredict which can classify a novel drug according to the 
Grandi 2010 APD as suggested in the paper. For this reason the Grandi APDs have been
added to the drug_data.dat file.

It should be compiled using
scons cl=1 exe=1 projects/CardiovascRes11/apps