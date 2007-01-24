#!/usr/bin/env python

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
