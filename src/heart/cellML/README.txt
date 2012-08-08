This folder contains files taken from the CellML repository (as of August 2010).

The -conf.xml files provide any extra information Chaste requires to translate them into C++.
(We add metadata modifiers and extend lookup tables as some of the pacing protocols push the 
models outside their normal ranges.)

The Shannon 2004 rabbit model is in the main Chaste repository.

The project uses Chaste to generate versions of the models which use Cvode as the solver, 
include partial-evaluation to optimise for speed, but make no simplifications such as look-up tables.

