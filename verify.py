#!/usr/bin/env python

import prime
n=prime.primorial(100)

f=open("p.txt","r")
for line in f:
  p, r = [long(i) for i in line.split()]
  cr = n%p
  if (r*cr)%p != 1:
    print "Fail", p, cr, r
