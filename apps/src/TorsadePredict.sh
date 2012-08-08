#!/bin/bash

# A wrapper script for Chaste that can figure out where it really lives on the
# filesystem, set the LD_LIBRARY_PATH to the correct subfolder and run Chaste.

#IFS=" \t\n"
#declare -x PATH=/bin:/usr/bin

script="$0"

# Have we been called via a symlink?
while [[ -L "$script" ]]; do script=$(readlink -n "$script"); done

# Figure out the folder $script is in
script_path=$(2>/dev/null cd "${script%/*}" >&2; echo "`pwd -P`/${script##*/}")
script_dir=$(dirname "$script_path")

#Set the LD_LIBRARY_PATH and run Chaste
export LD_LIBRARY_PATH="$script_dir/libs"

#Inform the user where the output will appear
if [ -z "$CHASTE_TEST_OUTPUT" ]; then
  echo "\$CHASTE_TEST_OUTPUT is currently unset.  Your output will appear in ./testoutput"
else
  echo "\$CHASTE_TEST_OUTPUT is currently set to " $CHASTE_TEST_OUTPUT. 
fi

# This line actually run Chaste with the given arguments
"$script_dir/TorsadePredict" "$@"
