// Computes 100# mod p
// p is read from a shared buffer and the result written back to the buffer

#include "e_lib.h"

#define NUM_LEN 3
#define EL_BITS 15
#define EL_MASK ((1 << EL_BITS) - 1)
typedef short el;
typedef int el2;

typedef struct inout_data
{
  long long p;
  long long res;
} data_t;
data_t buf[256] SECTION("shared_dram");

typedef el num_t[NUM_LEN];
typedef el numl_t[NUM_LEN+1];

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
  numl_t t;
  numl_t r = {0};
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
  }
  set(res, r);
}

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

num_t a, b, c;

unsigned long long doprimorial(unsigned long long p64)
{
  num_t p, abar, bbar, res, mi, t;
  el invp = invpow2(p64, 1 << EL_BITS);
  sets(p, p64);
  mont_init(mi, p64);

  mont_step(t, a, mi, p, invp);
  set(abar, t);
  mont_step(t, b, mi, p, invp);
  set(bbar, t);
  mont_step(t, abar, bbar, p, invp);
  mont_step(res, t, c, p, invp);

  return (unsigned long long)res[2] << (EL_BITS*2) | 
         (unsigned long long)res[1] << (EL_BITS) |
         res[0];
}

int main(void)
{
  sets(a, 223092870ll);
  sets(b, 2756205443ll);
  sets(c, 907383479ll);

  for (int i = 0; i < 256; ++i)
  {
    buf[i].res = doprimorial(buf[i].p);
  }

  __asm__ __volatile__ ("idle");
  return 0;
}

