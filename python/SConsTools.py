"""Useful functions for use by the build system."""

import os

def IsTemplateCpp(filepath): 
    """Check if the .cpp file defines a templated class. 
     
    We assume this is the case if one of the first 2 lines starts '#ifndef '. 
    """ 
    fp = open(filepath) 
    line = fp.next() 
    template = line.startswith('#ifndef ') 
    if not template: 
        line = fp.next() 
        template = line.startswith('#ifndef ') 
    fp.close() 
    return template

def FindSourceFiles(rootDir, ignoreDirs=[]):
    """Look for source files under rootDir.
    
    Returns 2 lists: the first of source (.cpp) files, and the second
    of the directories in which they may be found.
    
    Optionally, specify ignoreDirs to not search within particular
    folder names.
    """
    source_files = []
    source_dirs = []
    ignoreDirs.append('.svn')
    for dirpath, dirnames, filenames in os.walk(rootDir):
        for dirname in dirnames[:]:
            if dirname in ignoreDirs:
                dirnames.remove(dirname)
            else:
                source_dirs.append(os.path.join(dirpath, dirname))
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            if filename[-4:] == '.cpp' and not IsTemplateCpp(filepath):
                source_files.append(filepath)
    return source_files, source_dirs
