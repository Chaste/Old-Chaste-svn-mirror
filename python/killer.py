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

# Kill off any processes run by our user in this directory
# (except us, of course)

import os
import signal
import sys

sim = '-s' in sys.argv
if '-d' in sys.argv:
  i = sys.argv.index('-d')
  kill_dir = os.path.realpath(sys.argv[i+1])
else:
  kill_dir = os.path.realpath(os.getcwd())

print "Killing processes owned by", os.getpid(), "in", kill_dir

for pid in os.listdir('/proc/'):
  if pid[0] in "0123456789":
    try:
      s = os.stat('/proc/' + pid)
      if s.st_uid == os.getuid():
        cwd = os.path.realpath('/proc/' + pid + '/cwd')
        if cwd == kill_dir and int(pid) != os.getpid():
          print pid, "is running from our dir as"
          os.system('cat /proc/' + pid + '/cmdline')
          if not sim:
            os.kill(int(pid), signal.SIGTERM)
            print " ** SENT SIGTERM  mwa ha ha"
    except OSError:
      # We can't read all our processes; that's ok
      pass
