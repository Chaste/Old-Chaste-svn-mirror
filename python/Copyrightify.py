#!/usr/bin/env python

"""Copyright (C) University of Oxford, 2008

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



# Check, apply or modify the copyright notice
import os, sys
exts = ['.cpp', '.hpp', '.py']
dir_ignores = ['build', 'cxxtest', 'testoutput', 'doc', 'anim', 'projects']
#exclusions = ['triangle/triangle.cpp']

apply_update =  '-update' in sys.argv
apply_new = '-new' in sys.argv

chaste_dir = '.'
if '-dir' in sys.argv: 
    i = sys.argv.index('-dir')
    chaste_dir = os.path.realpath(sys.argv[i+1])


depricated_notice="""Copyright (C) University of Oxford, 2008

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
along with Chaste.  If not, see <http://www.gnu.org/licenses/>."""


current_notice="""Copyright (C) University of Oxford, 2008

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

#py_depricated_notice=''
#for line in depricated_notice.splitlines():#
#	py_depricated_notice+=''.join(['# ',line,'\n'])
#py_current_notice=''
#for line in current_notice.splitlines():
#	py_current_notice+=''.join(['# ',line,'\n'])

py_current_notice='\"\"\"'+current_notice+'\"\"\"\n'
cpp_current_notice='/*\n\n'+current_notice+'\n*/'
cpp_depricated_notice='/*\n'+depricated_notice+'\n*/'


pycml_notice="// Processed by pycml - CellML Tools in Python"
xsd_notice="// Copyright (C) 2005-2007 Code Synthesis Tools CC"
triangle_notice="""/*  Copyright 1993, 1995, 1997, 1998, 2002, 2005                             */
/*  Jonathan Richard Shewchuk                                                */"""
recursive_notice="___This file here___"

def CheckForCopyrightNotice(findStr, fileIn):
    fileIn.seek(0)
    file_text=fileIn.read()
    return (file_text.find(findStr) >= 0)
    
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
    if (fileName[-15:]=='Copyrightify.py'):
    	#Can't really check this one, since it knows all the licences
    	return True
    valid_notice=False
    if (CheckForCopyrightNotice(cpp_current_notice, file_in) or CheckForCopyrightNotice(py_current_notice, file_in)):
        #print 'Found current notice in '+file_name
        valid_notice=True
    if (CheckForCopyrightNotice(recursive_notice, file_in) or CheckForCopyrightNotice(pycml_notice, file_in) or CheckForCopyrightNotice(xsd_notice, file_in) or CheckForCopyrightNotice(triangle_notice, file_in)):
        #print 'Found 3rd party notice in '+file_name
        if (valid_notice):
            print "Multiple notices on"+file_name
            return False
        else:
            return True
    if (valid_notice):
        return True
    if (CheckForCopyrightNotice(cpp_depricated_notice, file_in)):
        print 'Found depricated copyright notice for',fileName
        if (apply_update):
           	ReplaceStringInFile(cpp_depricated_notice, cpp_current_notice, fileName)
        	return True
        else:
            print 'Fix this by running with -update argument'
            return False
    
    print 'Found no copyright notice for',fileName
    if (apply_new):
        if (fileName[-3:] == '.py'):
          	print 'Not implemented'
          	return False
        else:
	        HeadAppendStringInFile(cpp_current_notice, fileName)
        return True
    else:
        print 'Fix this by running with -new argument'
        return False

#os.chdir(chaste_dir)
num_no_copyrights=0
num_copyrights=0
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
            if (InspectFile(file_name) == False):
                num_no_copyrights+=1
            else:
                num_copyrights+=1

# Let the test summary script know
if num_no_copyrights > 0:
    print
    print "The next line is for the benefit of the test summary scripts."
    print "Failed",num_no_copyrights,"of",num_no_copyrights+num_copyrights,"tests"

    # Return a non-zero exit code if orphans were found
    import sys
    sys.exit(num_no_copyrights)
else:
    print "Copyright test passed ok."
print "Copyright test run over "+chaste_dir