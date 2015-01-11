#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <e-hal.h>

#define ENTRIES_PER_CORE 32768
#define NUM_CORES 16
#define ENTRIES (ENTRIES_PER_CORE*NUM_CORES)
typedef struct shared_data
{
  long long p[ENTRIES];
  long long res[ENTRIES];
} data_t;

#define _BufSize   (sizeof(buf))
#define _BufOffset (0x01000000)

data_t buf, buf2;

#define NUM_LOW_PRIMES 26000
unsigned int low_primes[NUM_LOW_PRIMES];
long long lastp = 0xffffffff;

//#define SIEVE_SIZE (256u*1024u*1024u*8u)
#define SIEVE_SIZE (256u*100u*8u)
unsigned int* mainsieve;

void genlowprimes()
{
  // Do something simple.
  unsigned int p, s, i;
  low_primes[0] = 3;
  low_primes[1] = 5;
  p = 7;
  s = 3;
  i = 2;
  while (i < sizeof(low_primes) / sizeof(int))
  {
    int j = 0;
    for (; low_primes[j] <= s; ++j)
    {
      if (p%low_primes[j] == 0)
        break;
    }
    if (low_primes[j] > s)
      low_primes[i++] = p;
    p += 2;
    if (s*s < p) ++s;
  }
}

static unsigned int genprimsieve[256*1024];
static unsigned int pattern357[105];
static unsigned genprimsievej;
static unsigned int genprimoffsets[NUM_LOW_PRIMES];
static long long base;

void fillsieve();

void initgenprimes()
{
  for (int i = 0; i < 3; ++i)
  {
    unsigned offset = 0;
    while (offset < sizeof(pattern357) * 8)
    {
      pattern357[offset >> 5] |= 1<<(offset&0x1f); 
      offset += low_primes[i];
    }
  }

  base = lastp+2;
  int i, j;
  unsigned max_prime=sqrt(base + sizeof(genprimsieve)*16)+1;
  for (i = 3; i < NUM_LOW_PRIMES; ++i)
  {
    unsigned offset = base % low_primes[i];
    if (offset) offset = low_primes[i] - offset;
    if (offset & 1) offset += low_primes[i];
    offset >>= 1;
    genprimoffsets[i] = offset;
  }
  base -= sizeof(genprimsieve)*16;
  fillsieve();
} 

void fillsieve()
{
  base += sizeof(genprimsieve)*16;
  memset(genprimsieve, 0, sizeof(genprimsieve));
  int i;
  unsigned max_prime=sqrt(base + sizeof(genprimsieve)*16)+1;
  
  unsigned offset = base % 105;
  if (offset != 0) {
    if (offset & 1) offset += 105;
    offset >>= 1;
    offset += 105*(offset * 7);
    offset >>= 5;
    offset %= 105;
  }
  for (i = 0; i < sizeof(genprimsieve) / sizeof(int); ++i)
  {
    genprimsieve[i] = pattern357[offset];
    if (++offset == 105) offset = 0;
  }
  for (i = 3; low_primes[i] < max_prime; ++i)
  {
    unsigned offset = genprimoffsets[i];
    while (offset < sizeof(genprimsieve) * 8)
    {
      genprimsieve[offset >> 5] |= 1<<(offset&0x1f); 
      offset += low_primes[i];
    }
    genprimoffsets[i] = offset - (sizeof(genprimsieve) * 8);
  }
  for (; i < NUM_LOW_PRIMES; ++i)
  {
    unsigned offset = genprimoffsets[i];
    unsigned mult = (sizeof(genprimsieve) * 8) / low_primes[i];
    offset += mult * low_primes[i];
    if (offset < sizeof(genprimsieve) * 8)
      offset += low_primes[i];
    genprimoffsets[i] = offset - (sizeof(genprimsieve) * 8);
  }
  genprimsievej = 0;
}

// Fill buf with the next primes
void genprimes(data_t* bufptr)
{
  int i, j;

  for (i = 0; i < ENTRIES; ++genprimsievej)
  {
    if (genprimsievej == sizeof(genprimsieve) * 8)
      fillsieve();
    unsigned j = genprimsievej;
    if ((genprimsieve[j>>5] & (1 << (j&0x1f))) == 0)
    {
      bufptr->p[i++] = base + (j << 1);
      //if ((bufptr->p[i-1] % 5) == 0)
      //  fprintf(stderr, "Sieve broken! %lld\n", bufptr->p[i-1]);
    }
  }
  lastp = bufptr->p[i-1];
}

// return t such that at = 1 mod m
unsigned inverse_short(unsigned a, unsigned m)
{
  int t = 0, newt = 1, r = m, newr = a;
  while (newr != 0)
  {
    int q = r / newr;
    int x = t - q * newt;
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

// Assumes m is odd.
unsigned inverse(unsigned a, unsigned b)
{
  if (b < 0x80000000)
    return inverse_short(a, b);

  unsigned alpha, beta;
  long long u, v, s, t;
  u = 1; v = 0; s = 0; t = 1;
  alpha = a; beta = b;

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

void sievep(unsigned p)
{
  unsigned n = (223092870ll * 2756205443ll) % p;
  n = (n * 907383479ll) % p;
  n = (n * 4132280413ll) % p;
  unsigned ni = inverse(n, p);
//  if (((unsigned long long)ni * n) % p != 1)
//    fprintf(stderr, "Bad inverse %u %u %u\n", ni, n, p);
  unsigned r = (ni == 0) ? 0 : (p - ni);
  for (; r < SIEVE_SIZE; r+=p)
    mainsieve[r>>5] |= 1<<(r&0x1f);  
}

void initsieve()
{
  unsigned i, j;
  for (i = 0; i < sizeof(low_primes)/sizeof(int); ++i)
  {
    if (low_primes[i] < 100) continue;
    sievep(low_primes[i]);
  }
  return;
  
#define INIT_SIEVE_MB 32
  unsigned base = low_primes[i-1]+2;
  unsigned int* sieve = malloc(INIT_SIEVE_MB*1024*1024);
  
  while (1)
  {
    memset(sieve, 0, INIT_SIEVE_MB*1024*1024);
    unsigned max_prime=sqrt(base + INIT_SIEVE_MB*1024*1024*16.0)+1;
    fprintf(stderr, "...\r");
    for (i = 0; low_primes[i] < max_prime; ++i)
    {
      unsigned offset = base % low_primes[i];
      if (offset) offset = low_primes[i] - offset;
      if (offset & 1) offset += low_primes[i];
      offset >>= 1;
      while (offset < INIT_SIEVE_MB*1024*1024 * 8)
      {
        sieve[offset >> 5] |= 1<<(offset&0x1f); 
        offset += low_primes[i];
      }
    }
    for (i = 0; i < INIT_SIEVE_MB*1024*1024*8 && (base+(i<<1)) >= base; ++i)
    {
      if ((sieve[i>>5] & (1 << (i&0x1f))) == 0)
      {
        sievep(base + (i << 1));
      }
      if ((i&0xffffff) == 0) 
      {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        fprintf(stderr, "%u %ld.%ld\n", base+(i<<1), ts.tv_sec, ts.tv_nsec/1000000);
      }
    }
    if (base+(i<<1) < base)
      break;
    base += i << 1;
    fprintf(stderr, "%u\n", base);
  }
  free(sieve);
}

void applysieve(data_t* bufptr)
{
  for (int i = 0; i < ENTRIES; ++i)
  {
    if (bufptr->res[i] == 0)
    {
      fprintf(stderr, "Bad: res=0 for %lld\n", bufptr->p[i]);
    }
    else
    {
      unsigned long long r = bufptr->p[i] - bufptr->res[i];
      if (r < SIEVE_SIZE)
      {
        if ((mainsieve[r>>5] & (1<<(r&0x1f))) == 0) printf("%lld %lld\n", r, bufptr->p[i]);
        mainsieve[r >> 5] |= 1<<(r&0x1f);
      }
    }
  }
}

int main(int argc, char*argv[])
{
  long long time = 0;
  int max_loops = 128;
  unsigned row, col, coreid, first;
  e_platform_t platform;
  e_epiphany_t dev;
  e_mem_t emem;
  data_t *curbuf = &buf, *nextbuf = &buf2;
  struct timespec sleeptime;
  sleeptime.tv_sec = 0;
  sleeptime.tv_nsec = 100000;
  
  // Generate low primes for sieving
  printf("Gen Low Primes\n");
  genlowprimes();

  // Allocate and init main sieve
  printf("Init main sieve\n");
  mainsieve = malloc(SIEVE_SIZE >> 3);
  memset(mainsieve, 0, SIEVE_SIZE >> 3);
  initsieve();

  // Start up epiphany
  printf("Init epiphany\n");
  e_init(NULL);
  e_reset_system();
  e_get_platform_info(&platform);

  // Allocate a buffer in shared external memory
  // for message passing from eCore to host.
  e_alloc(&emem, _BufOffset, _BufSize);

  // Open a workgroup
  e_open(&dev, 0, 0, platform.rows, platform.cols);

  // Load the device program onto all the eCores
  e_load_group("e_primorial.srec", &dev, 0, 0, platform.rows, platform.cols, E_FALSE);

  // Generate first batch of primes to work on
  printf("Init gen primes\n");
  initgenprimes();
  genprimes(nextbuf);
  first = 1;

  // Main loop
  for (int loops = 0; loops < max_loops; ++loops)
  {
    unsigned sleeps = 0;
    unsigned long long start, end; 
    struct timespec ts;

    // Copy in next batch of primes for epiphany and start the cores
    e_write(&emem, 0, 0, 0x0, nextbuf->p, sizeof(buf.p));

    clock_gettime(CLOCK_MONOTONIC, &ts);
    start = ts.tv_sec * 1000000000ll + ts.tv_nsec;

    for (int i = 0; i < 16; ++i)
    {
      e_start(&dev, i>>2, i&3);
    }

    // Apply last results (if not first time)
    if (!first)
      applysieve(curbuf);
    else
      first = 0;

    // Swap buffers
    {
      data_t* tmp = nextbuf;
      nextbuf = curbuf;
      curbuf = tmp;
    }

    // Generate next batch of primes
    genprimes(nextbuf);

    // Wait for epiphany cores.
    for (int i = 0; i < 16; ++i)
    {
      while (1)
      {
        unsigned status;
        e_read(&dev, i>>2, i&3, E_REG_STATUS, &status, sizeof(status));
        if ((status & 1) == 0)
          break;
        nanosleep(&sleeptime, NULL);
        ++sleeps;
        if (sleeps == 10000) fprintf(stderr, "Core %d stuck\n", i);
      }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts);
    end = ts.tv_sec * 1000000000ll + ts.tv_nsec;
    end -= start;
    time += end;

    e_read(&emem, 0, 0, sizeof(buf.p), curbuf->res, sizeof(buf.res));

    {
      fprintf(stderr, "%lld %lld.%04lld %d\n", curbuf->p[ENTRIES-1], end / 1000000000, (end / 100000) % 10000, sleeps);
    }
  
#if 0
    for (int i = 0; i < 16; ++i)
      for (int j = 0; j < _BufEntriesPerCore; j+=255)
        printf("%lld %lld\n", buf2[i*_BufEntriesPerCore+j].p, buf2[i*_BufEntriesPerCore+j].res);
#endif
  }

  applysieve(curbuf);

#if 0
  for (unsigned i = 0; i < SIEVE_SIZE; ++i)
  {
    if ((mainsieve[i>>5] & (1<<(i&0x1f))) == 0)
      printf("%u\n", i);
  }
#endif

  printf("Total time: %lld.%03lld\n", time / 1000000000, (time / 1000000) % 1000);
  time /= max_loops;
  printf("Time per loop: %lld.%04lld\n", time / 1000000000, (time / 100000) % 10000);

  e_close(&dev);
  e_free(&emem);
  e_finalize();

  return 0;
}
