#!/usr/bin/env python
import prime

primes = prime.sundaram3(100)

a=1
for p in primes:
  if a*p < (1<<32):
    a = a*p
  else:
    print a
    a = p
print a
