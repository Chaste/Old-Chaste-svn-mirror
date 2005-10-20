#!/bin/bash

# Run a `nightly' build on the integration machine.

# Pass a build type as first argument to build that, otherwise will
# build a default selection.

# The integration machine
machine='userpc60.comlab.ox.ac.uk'
# Base path for where test results will appear
remote_base='testresults/'
# Base path where test results should be put
local_base='/var/www/chaste_test_data/nightly/'

# Latest revision
rev=`svnlook youngest /var/svn/chaste`

# Create output dir
if [ ! -d $local_base$rev ]; then
	mkdir $local_base$rev
fi

do_build ()
{
	# Run the build
	ssh bob@$machine ./builder $rev $1 2>&1 >/var/www/html/out/nightly_$rev.$machine.$1
	# Remove any old results from this revision, machine & build type
	/bin/rm -rf $local_base$rev/$machine.$1
	# Copy results
	scp -r bob@$machine:$remote_base$machine.$1 $local_base$rev/$machine.$1
	# Add a link to the build log
	ln -s /var/www/html/out/nightly_$rev.$machine.$1 $local_base$rev/$machine.$1/build.log
	# And remove from the integration machine
	ssh bob@$machine /bin/rm -r $remote_base$machine.$1
}

if [ -z "$1" ]; then
	# Nightly test pack (don't include default tests)
	do_build default_onlytests_Nightly
	# Memory tests
	do_build MemoryTesting
	do_build ParallelMemoryTesting
	do_build Parallel_Continuous
else
	do_build $1
fi
