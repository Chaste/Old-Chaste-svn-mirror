Using CellML with (Cardiac) Chaste
==================================

There is a companion tool to Chaste, named PyCml, which can generate
Chaste-compatible C++ code from cardiac ionic cell models described in
CellML.  This allows Chaste to make use of any such model.  It also
applies some optimisations to improve the speed of simulations.

Since PyCml generates C++ source code, the source code release
of Chaste is required to make use of CellML models.

Relevant information can also be found on the developers' wiki at
https://chaste.comlab.ox.ac.uk/cgi-bin/trac.cgi/wiki/CodeGenerationFromCellML
For a guest login, use the username "anonymous", and your email address as
the password.


Use of CellML in the cardiac executable
=======================================

As of release 2.0 of Chaste, the cardiac executable has gained the ability to
automatically load cell models encoded as CellML files at run-time, rather than
needing them to be incorporated within Chaste when it is compiled.  In order to
take advantage of this, you need (at present) to have built the executable from
source yourself, as it uses your Chaste source tree to convert the CellML file
into runnable code.

This change means that the definition of ionic models within the parameters file
has changed.  To specify one of the models included within Chaste, you now need to
wrap it in a <Hardcoded> element, e.g.

    <Hardcoded>FaberRudy2000</Hardcoded>

Specify a dynamically loaded model using the <Dynamic> element instead, e.g.

    <Dynamic>
        <Path relative_to="chaste_source_root">heart/dynamic/libDynamicallyLoadableLr91.so</Path>
    </Dynamic>

The path may either point to a pre-compiled shared library (as in the example),
or a CellML file.  If the latter, a shared library will be created in the same
folder as the CellML file, and loaded by the executable at run-time. See
heart/test/data/xml/ChasteParametersFullFormat.xml for a full example parameters
file.


Installing PyCml
================

PyCml itself is a pure Python tool distributed with the Chaste source
code (in folder python/pycml), but it does have two prerequisites:
Amara and RNV.

If you install the Ubuntu Chaste package, this should pull in everything
needed.

Instructions for manually installing the dependencies can be found in
INSTALLATION.txt, and also online at
https://chaste.comlab.ox.ac.uk/cgi-bin/trac.cgi/wiki/InstallPyCml


Using PyCml
===========

The process of code generation is mostly automatic, although some of
the options require some human intervention.  There is a helper script
included in the Chaste source release (python/ConvertCellModel.py) to
simplify the process of using PyCml.

Using the helper script
-----------------------

The script python/ConvertCellModel.py should ease the process of using PyCml.
Run the script supplying CellML files as arguments. Extra options for PyCml
can also be supplied. It will generate both normal and 'fully' optimised
versions of models (i.e. with both partial evaluation and lookup tables done).

Examples:
  python/ConvertCellModel.py --assume-valid heart/src/odes/cellml/luo_rudy_1991.cellml
  python/ConvertCellModel.py heart/src/odes/cellml/AnotherCellModel.cellml --no-member-vars

Some options are specified implicitly. At the time of writing these are:
 --conf=config.xml --use-chaste-stimulus --convert-interfaces --Wu --row-lookup-method

By default the helper script will use the PyCml shipped with Chaste. If,
however, you have the environment variable PYCML_DIR set, the script will
assume this points to an installation of PyCml that you want to use instead. 

The PyCml configuration file
----------------------------

Various parts of PyCml's functionality will not work properly without
a configuration file being specified.  Notably this is needed to
indicate which variables in the model represent the transmembrane
potential, stimulus current, and other transmembrane ionic currents.
The configuration file also specifies which variable(s) are used to
index lookup tables, and the ranges over which these vary.

The configuration file is XML, and contains both general and
model-specific settings.  A sample file is provided with PyCml.  This
may be modified to add new model-specific configuration, or change the
global settings.

Further information on the configuration file is given in the
documentation for ConfigurationStore.read_configuration_file in
translate.py.  The key settings are as follows.

Global options are specified in a global element, for example:

<global>
  <lookup_tables>
    <lookup_table>
      <var>membrane,V</var> <!-- The index variable for the table -->
      <min>-100.0001</min>
      <max>49.9999</max>
      <step>0.01</step> <!-- The table step size -->
    </lookup_table>
  </lookup_tables>
  <currents>
    <stimulus>membrane,i_Stim</stimulus>
    <ionic_match>membrane,i_.*</ionic_match> <!-- regexp on var names -->
    <!--
      Note that the stimulus current is automatically
      excluded from being an ionic current.
      Also note that there are implicit ^ and $ around the regexp.
      -->
  </currents>
  <transmembrane_potential>membrane,V</transmembrane_potential>
</global>

Note that variable names are given in full form,
i.e. "component,variable".  The component given should be that from
which the variable is exported.  Ionic currents are specified using a
regular expression to match variable names.

Model-specific settings are given in for_model or for_models elements.
They can contain the same options as in the global element.  Models
can be specified by name (in the case of for_model) or by id (in both
cases).  For example:

<for_model id="mahajan_shiferaw_2008_version01">
  ...
</for_model>
<for_model name="fox_model_2002">
  ...
</for_model>
<for_models>
  <ids>
    <id>ten_tusscher_model_2004_endo</id>
    <id>ten_tusscher_model_2004_M</id>
    <id>ten_tusscher_model_2004_epi</id>
  </ids>
  ...
</for_models>

Invoking PyCml directly
-----------------------

An important point to note is that, at present, the tools must be run with the
current working directory set to their location, i.e. cd /path/to/pycml before
running the commands below.

Another general point is that, with any of the transformations, you may need to
include the --Wu flag if the model does not pass units validation. This demotes
the error messages to warnings, allowing transformation to proceed (assuming no
other errors were found). 

> Generating a standard AbstractCardiacCell subclass

This one is easy, and the most automatic.

./translate.py --conf=config.xml --use-chaste-stimulus --convert-interfaces <cellml_file>

> Generating optimised AbstractCardiacCell subclasses

There is a helper script shipped with PyCml to generate each
combination of optimisations.

./do_trans.py --use-chaste-stimulus --convert-interfaces <cellml_file>

If you just want full optimisation, use

./translate.py --conf=config.xml --use-chaste-stimulus --convert-interfaces -a -p -l <cellml_file>

> Generating an AbstractBackwardEulerCardiacCell subclass

This is a multi-stage process, as the analysis uses Maple to perform
symbolic differentiation.  You will thus need access to an
installation of Maple to use this functionality.

./translate.py -J --omit-constants --conf=config.xml <cellml_file.cellml>
maple -i <cellml_file.mpl> > <cellml_file.out>
./translate.py -j <cellml_file.out> --conf=config.xml --use-chaste-stimulus --convert-interfaces <cellml_file.cellml>

> Generating an optimised AbstractBackwardEulerCardiacCell subclass

The following commands will apply all optimisations.  Note that the
lookup tables analysis & partial evaluation are not yet applied to
some portions of the model, notably the Jacobian calculation (although
Maple does perform much PE for us) and linear updates in
ComputeExceptVoltage.  Hence the optimisations don't make much
difference yet.

./translate.py -J --conf=config.xml <cellml_file.cellml>
maple -i <cellml_file.mpl> > <cellml_file.out>
./translate.py -j <cellml_file.out> --conf=config.xml --use-chaste-stimulus --convert-interfaces -a -p -l <cellml_file.cellml>
