// Computes 100#^-1 mod p
// p is read from a shared buffer and the result written back to the buffer

#include "e_lib.h"

#define NUM_LEN 3
#define EL_BITS 15
#define EL_MASK ((1 << EL_BITS) - 1)
#define ENTRIES_PER_CORE 32768
#define NUM_CORES 16
#define ENTRIES (ENTRIES_PER_CORE*NUM_CORES)

typedef int el; // Actually use int not short to avoid conversions
                // when going in and out of registers.
typedef int el2;


typedef struct shared_data
{
  long long p[ENTRIES];
  long long res[ENTRIES];

} data_t;
data_t* buf = (data_t*)0x8f000000;

typedef el num_t[NUM_LEN];
typedef el numl_t[NUM_LEN+1];

static void zero(el* v)
{
  for (int i = 0; i < NUM_LEN; ++i)
    v[i] = 0;
}

static void sets(el* a, long long b)
{
  for (int i = 0; i < NUM_LEN; ++i)
  {
    a[i] = (el)(b & EL_MASK);
    b >>= EL_BITS;
  }
}

static void set(el* a, el* b)
{
  for (int i = 0; i < NUM_LEN; ++i)
    a[i] = b[i];
}

static void addl(el* res, el* a, el* b)
{
  el c = 0;
  for (int i = 0; i < NUM_LEN+1; ++i)
  {
    el2 t = a[i] + b[i] + c;
    c = (el)(t >> EL_BITS);
    res[i] = (el)(t & EL_MASK);
  }
} 

static void mult(el* res, el* a, el b)
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

// res = (a*b) + c
static void mult_add(el* res, el* a, el b, el* c)
{
  el d = 0;
  for (int i = 0; i < NUM_LEN; ++i)
  {
    el2 t = a[i] * b + c[i] + d;
    d = (el)(t >> EL_BITS);
    res[i] = (el)(t & EL_MASK);
  }
  res[NUM_LEN] = c[NUM_LEN] + d;
}

static void mont_step(el* res, el* a, el* b, el* p, el invp)
{
  numl_t r = {0};
  for (int i = 0; i < NUM_LEN; ++i)
  {
    mult_add(r, a, b[i], r);
    el2 q = ((1 << EL_BITS) - r[0]) * invp;
    mult_add(r, p, (el)(q & EL_MASK), r);
    for (int j = 1; j <= NUM_LEN; ++j)
    {
      r[j-1] = r[j];
    }
    r[NUM_LEN] = 0;
  }
  set(res, r);
}

static void mont_init(el* res, unsigned long long p)
{
  unsigned rhigh = 1;
  unsigned phigh = p >> 32;
  unsigned long long r;
  int i;
  for (i = 32; rhigh < phigh; ++i)
    rhigh <<= 1;
  r = (unsigned long long)rhigh << 32;
  for (; i < 90; ++i)
  {
    if (r > p) r -= p;
    r <<= 1;
  }
  if (r > p) r -= p;
  sets(res, r);
}

static unsigned long long inverse(unsigned long long a, unsigned long long b)
{
  long long alpha, beta;
  long long u, v, s, t;
  u = 1; v = 0; s = 0; t = 1;
  alpha = a; beta = b;

  if (a == 0)
    return 0;

  // Keep a = u * alpha + v * beta
  while ((a&1) == 0)
  {
    a >>= 1;
    if ((u|v) & 1)
    {
      u = (u + beta) >> 1;
      v = (v - alpha) >> 1;
    }
    else
    {
      u >>= 1;
      v >>= 1;
    }
  }
  while (a!=b)
  {
    if ((b&1)==0)
    {
      b >>= 1;
      if ((s|t) & 1)
      {
        s = (s + beta) >> 1;
        t = (t - alpha) >> 1;
      }
      else
      {
        s >>= 1;
        t >>= 1;
      }
    }
    else if (b < a)
    {
      long long tmp;
      tmp = a;
      a = b;
      b = tmp;
      tmp = u;
      u = s;
      s = tmp;
      tmp = v;
      v = t;
      t = tmp;
    }
    else
    {
      b = b - a;
      s = s - u;
      t = t - v;
    }
  }
  if (a > 1) return 0;
  while (s < 0) s += beta;
  while (s >= beta) s -= beta;
  return s;
}

static unsigned invpow2(unsigned a, unsigned m)
{
  // m is a power of 2.
  m >>= 1;
  unsigned u = 1, v = 0;
  unsigned alpha = m, beta = a;

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

#define NUM_PRIM 4
num_t prim[NUM_PRIM];

unsigned long long doprimorial(unsigned long long p64)
{
  num_t p, primbar, res, mi, t, u, one;
  el invp = invpow2(p64 & EL_MASK, 1 << EL_BITS);
  sets(p, p64);
  mont_init(mi, p64);

  // Compute t = n# (mont space)
  mont_step(t, prim[0], mi, p, invp);
  for (int i = 1; i < NUM_PRIM-1; ++i)
  {
    mont_step(primbar, prim[i], mi, p, invp);
    mont_step(u, primbar, t, p, invp);
    set(t, u);
  }
  
  mont_step(res, t, prim[NUM_PRIM-1], p, invp);

  unsigned long long r64 =
         (unsigned long long)res[2] << (EL_BITS*2) | 
         (unsigned long long)res[1] << (EL_BITS) |
         res[0];
  while (r64 >= p64) r64 -= p64;

  r64 = inverse(r64, p64);

  return r64;
}

void null_isr(int);

int main(void)
{
  e_coreid_t coreid;
  coreid = e_get_coreid();
  unsigned int row, col, core;
  e_coords_from_coreid(coreid, &row, &col);
  core = row * 4 + col;

  sets(prim[0], 223092870ll);
  sets(prim[1], 2756205443ll);
  sets(prim[2], 907383479ll);
  sets(prim[3], 4132280413ll);

  // Must set up a null interrupt routine to allow host to wake us from idle
  e_irq_attach(E_SYNC, null_isr);
  e_irq_mask(E_SYNC, E_FALSE);
  e_irq_global_mask(E_FALSE);

  unsigned long long* primeptr = &buf->p[core*ENTRIES_PER_CORE];
  unsigned long long* resptr = &buf->res[core*ENTRIES_PER_CORE];

  while (1)
  {
    // Always work to do immediately
    for (int i = 0; i < ENTRIES_PER_CORE; ++i)
    {
      resptr[i] = doprimorial(primeptr[i]);
    }

    // Wait for more work
    __asm__ __volatile__ ("idle");
  }

  return 0;
}

void __attribute__((interrupt)) null_isr(int x) { return; }
