def sundaram3(max_n):
    numbers = range(3, max_n+1, 2)
    half = (max_n)//2
    initial = 4

    for step in xrange(3, max_n+1, 2):
        for i in xrange(initial, half, step):
            numbers[i-1] = 0
        initial += 2*(step+1)

        if initial > half:
            return [2] + filter(None, numbers)

sieve = sundaram3(70000)

def nextprime(p):
  p+=2
  while True:
    for q in sieve:
      if p%q == 0:
        p+=2
        break
    else:
      return p

