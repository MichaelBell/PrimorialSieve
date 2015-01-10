#include <stdio.h>

#define NUM_LEN 3
#define EL_BITS 15
#define EL_MASK ((1 << EL_BITS) - 1)
typedef short el;
typedef int el2;

void zero(el* v)
{
  for (int i = 0; i < NUM_LEN; ++i)
    v[i] = 0;
}

void sets(el* a, long long b)
{
  for (int i = 0; i < NUM_LEN; ++i)
  {
    a[i] = (el)(b & EL_MASK);
    b >>= EL_BITS;
  }
}

void set(el* a, el* b)
{
  for (int i = 0; i < NUM_LEN; ++i)
    a[i] = b[i];
}

void addl(el* res, el* a, el* b)
{
  el c = 0;
  for (int i = 0; i < NUM_LEN+1; ++i)
  {
    el2 t = a[i] + b[i] + c;
    c = (el)(t >> EL_BITS);
    res[i] = (el)(t & EL_MASK);
  }
} 

void mult(el* res, el* a, el b)
{
  el c = 0;
  for (int i = 0; i < NUM_LEN; ++i)
  {
    el2 t = a[i] * b + c;
    c = (el)(t >> EL_BITS);
    res[i] = (el)(t & EL_MASK);
  }
  res[NUM_LEN] = c;
} 

void mont_step(el* res, el* a, el* b, el* p, el invp)
{
  el t[NUM_LEN+1];
  el r[NUM_LEN+1] = {0};
  for (int i = 0; i < NUM_LEN; ++i)
  {
    mult(t, a, b[i]);
    addl(r, r, t);
    el2 q = ((1 << EL_BITS) - r[0]) * invp;
    mult(t, p, (el)(q & EL_MASK));
    addl(r, r, t);
    for (int j = 1; j <= NUM_LEN; ++j)
    {
      r[j-1] = r[j];
    }
    r[NUM_LEN] = 0;
    printf("%d %d\n", r[1], r[0]);
  }
  set(res, r);
}

// Hack
unsigned long long mulmod64(unsigned long long a, unsigned long long b, unsigned long long p)
{
  unsigned long long acc = 0;
  for (int i = 45; i >= 0; --i)
  {
    acc <<= 1;
    acc = (acc > p) ? (acc - p) : acc;
    if ((b >> i) & 1)
    {
      acc += a;
      acc = (acc > p) ? (acc - p) : acc;
    }
  }
  return acc;
}

void mont_init(el* res, unsigned long long p)
{
  unsigned long long r = 1ll << 45;
  r %= p;
  r = mulmod64(r, r, p);
  sets(res, r);
}

// return t such that at = 1 mod m
unsigned long long inverse(unsigned long long a, unsigned long long m)
{
  unsigned long long t = 0, newt = 1, r = m, newr = a;
  while (newr != 0)
  {
    unsigned long long q = r / newr;
    unsigned long long x = t - q * newt;
    t = newt;
    newt = x;
    x = r - q * newr;
    r = newr;
    newr = x;
  }
  if (r > 1) return 0;
  if (t < 0) t += m;
  return t;
}

unsigned long long invpow2(unsigned long long a, unsigned long long m)
{
  // m is a power of 2.
  m >>= 1;
  unsigned long long u = 1, v = 0;
  unsigned long long alpha = m, beta = a;

  while (m > 0) 
  {
    m >>= 1;
    if ((u & 1) == 0)
    {
      u >>= 1;
      v >>= 1;
    }
    else
    {
      u = (u + beta) >> 1;
      v = (v >> 1) + alpha;
    }
  }

  return (alpha << 1) - v;
} 

void main()
{
  unsigned long long p64 = 7748603;
  unsigned long long a64 = 62342;
  unsigned long long b64 = 75558;
  el p[NUM_LEN], a[NUM_LEN], b[NUM_LEN], res[NUM_LEN], mi[NUM_LEN], t[NUM_LEN];
  //el invp = inverse(p64, 1 << EL_BITS);
  el invp = invpow2(p64, 1 << EL_BITS);
  sets(p, p64);
  mont_init(mi, p64);
  sets(a, a64);
  sets(b, b64);

  printf("%d %d %d %d\n", p[0], invp, mi[1], mi[0]);
  mont_step(t, a, mi, p, invp);
  set(a, t);
  printf("%d\n", a[0]);
  //mont_step(t, b, mi, p, invp);
  //set(b, t);
  //printf("%d\n", b[0]);
  mont_step(res, a, b, p, invp);
  //printf("%d\n", res[0]);
  //sets(a, 1);
  //mont_step(t, res, a, p, invp);

  printf("%d\n", (((el2)res[1])<<EL_BITS) + res[0]);
  printf("%lld\n", (a64*b64)%p64);
}

