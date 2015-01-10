#!/usr/bin/env python

from operator import mul

import prime
n=reduce(mul, prime.sundaram3(100), 1)

f=open("p.txt","r")
for line in f:
  p, r = [long(i) for i in line.split()]
  cr = n%p
  if r != cr:
    print "Fail", p, cr, r, r-p
