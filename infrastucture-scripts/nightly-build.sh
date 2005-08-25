#!/bin/bash

# Run a `nightly' build on the integration machine.

# Pass a build type as first argument to build that, otherwise will
# build a default selection.

# The integration machine
machine='userpc60.comlab.ox.ac.uk'
# Base path for where test results will appear
remote_base='eclipse/workspace/trunk/testoutput/'
# Base path where test results should be put
local_base='/var/www/chaste_test_data/nightly/'

# Latest revision
rev=`svnlook youngest /var/svn/chaste`

do_build ()
{
	# Run the build, discarding output
	ssh bob@$machine ./builder $rev $1 >/dev/null
	# Copy results
	scp -r bob@$machine:$path_base$machine.$1 $local_base$rev/$machine.$1
}

if [ -z "$1" ]; then
	# Memory tests
	do_build MemoryTesting
else
	do_build $1
fi
