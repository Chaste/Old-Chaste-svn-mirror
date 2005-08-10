#!/usr/bin/python

print "101    1     0    1"
for i in range(101):
  if(i == 0 or i == 100):
    print i, i/100.0, 1
  else:
    print i, i/100.0, 0
