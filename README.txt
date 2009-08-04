========= WELCOME TO CHASTE ===========

The files you have downloaded contain the source code for all the
Chaste functionalities.  Chaste makes use of a variety of external
libraries and packages that need to be installed on your machine.  The
file docs/INSTALLATION.txt provides a comprehensive guide on how to
install these external tools.

Chaste is distributed under the GNU Lesser General Public License.
The text of this licence is distributed in the file docs/Copying.txt.
Chaste uses various third party libraries which have their own
licences.  For details of these licences and the impact they may have
on your use of Chaste please see docs/Licences.html.

Chaste includes a complete test suite covering all the source
code. The easiest way to use existing source codes is to create a test
file which can call upon any of the source files, scons will build
this file for you and handle all of the dependencies and library
calls: e.g.
scons test_suite=projects/example/test/TestHello.hpp

We suggest you use the projects directory in this manner to store your
own source and test files if you do not wish to modify the chaste
source code.

For more information please refer to the chaste website at: 
http://web.comlab.ox.ac.uk/chaste/

Information on changes in this release can be found in the file
docs/ReleaseNotes.html.

Tutorial examples for this release are available at:
https://chaste.ediamond.ox.ac.uk/chaste/tutorials/release_1_1

Documentation generated from the code by Doxygen is available at:
https://chaste.ediamond.ox.ac.uk/chaste/docs/release_1_1
