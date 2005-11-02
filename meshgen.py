#!/usr/bin/env python

import os,sys

if len(sys.argv) < 2:
  sys.exit(1)

pref = sys.argv[1]
eles = 20
nodes = eles+1

f=file(pref+'.node', 'w')
f.write("%d\t1\t0\t1\n" % (nodes))
for i in range(0,nodes):
  b = 0
  if i == 0 or i == nodes: b=1
  f.write("%d %f %d\n" % (i, 1.0*i/eles, b))
f.close()

f=file(pref+'.ele', 'w')
f.write("%d\t2\t0\n" % (eles))
for i in range(0, eles):
  f.write("%d %d %d\n" % (i, i, i+1))
f.close()
