// Computes 100# mod p
// p is read from a shared buffer and the result written back to the buffer

#include "e_lib.h"

#define NUM_LEN 3
#define EL_BITS 15
#define EL_MASK ((1 << EL_BITS) - 1)
#define ENTRIES_PER_CORE 1024
typedef short el;
typedef int el2;

typedef struct inout_data
{
  long long p;
  long long res;
} data_t;
data_t buf[ENTRIES_PER_CORE*16] SECTION("shared_dram");

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

void mont_inverse(el* res, el* abar, unsigned long long p64, el* p, el invp)
{
  num_t base,t;
  unsigned long long exp = (p64 - 2) >> 1;
  set(res, abar);
  set(base, abar);
  while (exp > 0)
  {
    mont_step(t, base, base, p, invp);
    set(base, t);
    if (exp & 1)
    {  
      mont_step(t, base, res, p, invp);
      set(res, t);
    }
    exp >>= 1;
  }
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

#define NUM_PRIM 4
num_t prim[NUM_PRIM];

unsigned long long doprimorial(unsigned long long p64)
{
  num_t p, primbar, res, mi, t, u, one;
  el invp = invpow2(p64, 1 << EL_BITS);
  sets(p, p64);
  mont_init(mi, p64);

  // Compute t = n# (mont space)
  mont_step(t, prim[0], mi, p, invp);
  for (int i = 1; i < NUM_PRIM; ++i)
  {
    mont_step(primbar, prim[i], mi, p, invp);
    mont_step(u, primbar, t, p, invp);
    set(t, u);
  }
  
  // Compute u = (n#)^-1 (mont space)
  mont_inverse(u, t, p64, p, invp);

  // Get out of mont space.
  sets(one, 1);
  mont_step(res, u, one, p, invp);

  unsigned long long r64 =
         (unsigned long long)res[2] << (EL_BITS*2) | 
         (unsigned long long)res[1] << (EL_BITS) |
         res[0];
  while (r64 >= p64) r64 -= p64;
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

  data_t* bufptr = &buf[core*ENTRIES_PER_CORE];

  while (1)
  {
    // Always work to do immediately
    for (int i = 0; i < ENTRIES_PER_CORE; ++i)
    {
      bufptr[i].res = doprimorial(bufptr[i].p);
    }

    // Wait for more work
    __asm__ __volatile__ ("idle");
  }

  return 0;
}

void __attribute__((interrupt)) null_isr(int x) { return; }
