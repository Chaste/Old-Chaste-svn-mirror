#!/usr/bin/python

import os, sys
#Read the log file from the previous doxygen run
os.chdir('/root/')
logfile=open('doxygen-buildver.log', 'r')
lastBuiltVersion = int(logfile.readline())
logfile.close()
#Check to see if the version number has changed
buildVersion = int(os.popen('svnlook youngest /var/svn/chaste').readline())
#If so, do stuff
if(buildVersion > lastBuiltVersion or "-f" in sys.argv):
  print 'Repository has been updated, running doxygen...'
  sys.stdout.flush()
  #Check out repository
  os.system('svn co file:///var/svn/chaste/trunk/ tmp123/')
  os.chdir('tmp123/')
  #Run doxygen
  exitcode = os.system('( cat Doxyfile ; echo "PROJECT_NUMBER=Build:: '+str(buildVersion)+'" ) | doxygen - 2>&1 | tee ../doxygen-output.log')
  if exitcode:
    print "Doxygen returned non-zero exit code."
  if not os.path.exists('doxygen/html'):
    print "What happened to our output????"
  #Copy the files to /var/www/html/
  os.system('rm -rf /var/www/html/docs/')
  os.system('cp -r doxygen/html/ /var/www/html/docs/')
  #Write a log file saying the version of documentation that has been written
  os.chdir('../')
  os.system('rm -rf tmp123/')
  logfile=open('doxygen-buildver.log', 'w')
  logfile.write(str(buildVersion))
  logfile.close()
else:
  print 'No change to repository since previous execution'
