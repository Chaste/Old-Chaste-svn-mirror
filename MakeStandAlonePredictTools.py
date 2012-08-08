#!/usr/bin/env python


"""Copyright (C) University of Oxford, 2005-2010

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
"""

import os
import sys
import shutil
import re
import tarfile
import subprocess

# This script is to make a 32-bit or 64-bit binary release of the current Chaste

#Function to copy and strip of the xsd location
def copy_strip(in_file_name, out_file_name, strip_path):
  out_file = open(out_file_name, 'w')
  for line in open(in_file_name):
    out_file.write(line.replace(strip_path, ''))
  out_file.close()

# Here is the new release number:
release_number='0_5'

arch = os.uname()[-1] #Architecture is the last in the 5-tuple produced by uname
if ( arch == 'x86_64'):
 arch_bits='64'
else:
 arch_bits='32'

base_name = 'PredictTools_' + release_number + '_' + arch_bits + 'bit'
directory_name = base_name + '/'
tar_name = base_name + '.tgz'

print 'Making a standalone Gary''s Tools release in ' + directory_name
   
#Compile the code and bail out if necessary
print 'Compiling dynamically-linked executable'
if os.system('scons build=GccOpt static=0 chaste_libs=1 exe=1 compile_only=1 projects/GaryM/apps'):
  print 'General build failure.  You aren\'t ready to release!'
  sys.exit(1)

#Make the new directory for the packet and put the executable/wrapper in place
if os.path.exists(directory_name):
  shutil.rmtree(directory_name)
os.mkdir(directory_name)
shutil.copy('projects/GaryM/apps/src/ApPredict', directory_name)
shutil.copy('projects/GaryM/apps/src/TorsadePredict', directory_name)
shutil.copy('projects/GaryM/apps/src/ApPredict.sh', directory_name)
shutil.copy('projects/GaryM/apps/src/TorsadePredict.sh', directory_name)

print 'Making a subdirectory of library dependencies'
found_chaste_libraries=False
ldd=subprocess.Popen( ['ldd', directory_name+'ApPredict'], stdout=subprocess.PIPE ).communicate()[0]
os.mkdir(directory_name+'libs')
#Look for "LIB WHITESPACE => WHITESPACE PATH WHITESPACE"
lib_location = re.compile(r'(\S*)\s*=>\s*(\S*)\s*')
for lib_pair in lib_location.findall(ldd):
  if lib_pair[1][0]=='(' or lib_pair[1]=='not':
    print 'No library found for ', lib_pair[0]
  elif lib_pair[0]=='libc.so.6' or lib_pair[0]=='libpthread.so.0':
    print 'Ignoring library ', lib_pair[0], 'for compatibility'
  else:
    if lib_pair[0]=='libheart.so':
      found_chaste_libraries=True  
    shutil.copy(lib_pair[1], directory_name+'libs')

if found_chaste_libraries==False:
  print 'Could not find Chaste libraries (e.g. libheart.so).  Please set LD_LIBRARY_PATH.'
  sys.exit(1)

#Copy the xml files (with the hard coded location erased)
copy_strip('ChasteParameters.xml', directory_name+'ChasteParameters.xml', 'heart/src/io/')
copy_strip('heart/test/data/xml/ChasteParametersFullFormat.xml', directory_name+'ChasteParametersFullFormat.xml', '../../../src/io/')
copy_strip('heart/test/data/xml/ChasteParametersResumeSimulationFullFormat.xml', directory_name+'ChasteParametersResumeSimulationFullFormat.xml', '../../../src/io/')
#Copy all xsd files
os.system('cp heart/src/io/ChasteParameters*.xsd ' + directory_name)
os.system('cp projects/GaryM/test/drug_data.dat ' + directory_name)

# Following not needed for Gary's Tools
##Copy archive convert
#shutil.copy('archive_convert.sed', directory_name)
##Copy the standalone documentation to the highest level
#shutil.copy('docs/README_STANDALONE.txt', directory_name+'README.txt')
##Copy the rest of the documentation - without Subversion or editor roll-backs
#shutil.copytree('docs', directory_name+'docs')
#shutil.rmtree(directory_name+'docs/.svn')
#shutil.rmtree(directory_name+'docs/licences/.svn')
for dirpath, dirnames, filenames in os.walk(directory_name):
  for filename in filenames:
    if filename[-1] == '~':
      os.remove(dirpath+'/'+filename)

print 'Making tar file ' + tar_name
tar = tarfile.open(tar_name, "w:gz")
tar.add(directory_name)
tar.close()

#print 'Moving '+ tar_name +' to the correct place for downloads.  You will be prompted three times for a password'
#os.system('scp '+ tar_name + ' chaste@chaste.comlab.ox.ac.uk:')
#os.system('ssh -t chaste@chaste.comlab.ox.ac.uk "sudo mv '+ tar_name +' /var/www/html/chaste/downloads/"')
