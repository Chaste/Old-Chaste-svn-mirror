#!/usr/bin/env python

"""Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>."""

# Check, apply or modify the copyright notice
import os, sys
exts = ['.cpp', '.hpp']
dir_ignores = ['build', 'cxxtest', 'testoutput', 'doc', 'anim']

apply_update =  '-update' in sys.argv
apply_new = '-new' in sys.argv

chaste_dir = '.'
if '-dir' in sys.argv: 
    i = sys.argv.index('-dir')
    chaste_dir = os.path.realpath(sys.argv[i+1])


depricated_notice="""/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/"""


current_notice="""/*

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
"""


def CheckForCopyrightNotice(findStr, fileIn):
    fileIn.seek(0)
    file_text=fileIn.read()
    return (file_text.find(findStr) == 0)
    
def ReplaceStringInFile(findStr,repStr,filePath):
   "replaces all findStr by repStr in file filePath"
   tempName=filePath+'~'
   input = open(filePath)
   output = open(tempName,'w')

   s=input.read()
   output.write(s.replace(findStr,repStr))
   output.close()
   input.close()
   os.rename(tempName,filePath)
   print 'Notice: replaced depricated copyright notice in '+filePath

def HeadAppendStringInFile(appendString, filePath):
   "adds appendStr to the top of file filePath"
   tempName=filePath+'~'
   input = open(filePath)
   output = open(tempName,'w')

   s=input.read()
   output.write(appendString)
   output.write(s)
   output.close()
   input.close()
   os.rename(tempName,filePath)
   print 'Notice: applied copyright notice in '+filePath
    

   

def InspectFile(fileName):
    file_in = open(fileName)
    if (CheckForCopyrightNotice(current_notice, file_in)):
        #print 'Found current notice in '+file_name
        return True
    if (CheckForCopyrightNotice(depricated_notice, file_in)):
        print 'Found depricated copyright notice for',fileName
        if (apply_update):
            ReplaceStringInFile(depricated_notice, current_notice, fileName)
            return True
        else:
            print 'Fix this by running with -update argument'
            return False
    
    print 'Found no copyright notice for',fileName
    if (apply_new):
        HeadAppendStringInFile(current_notice, fileName)
        return True
    else:
        print 'Fix this by running with -new argument'
        return False

#os.chdir(chaste_dir)
num_no_copyrights=0
num_source_files=0
# for root, dirs, files in os.walk(chaste_dir):
for root, dirs, files in os.walk(chaste_dir):
    # Check for ignored dirs
    for dir in dir_ignores:
        if dir in dirs:
            dirs.remove(dir)
    # Check for source files
    for file in files:
        name, ext = os.path.splitext(file)
        if ext in exts:
            # Run astyle
            file_name = os.path.join(root, file)
            #ReplaceStringInFile('Chaste','Chaste2',os.path.join(root, file))
            num_source_files+=1
            if (InspectFile(file_name) == False):
                num_no_copyrights+=1

# Let the test summary script know
if num_no_copyrights > 0:
    print
    print "The next line is for the benefit of the test summary scripts."
    print "Failed",num_no_copyrights,"of",num_source_files,"tests"

    # Return a non-zero exit code if orphans were found
    import sys
    sys.exit(num_no_copyrights)
else:
    print "Copyright test passed ok."
print "Copyright test run over "+chaste_dir