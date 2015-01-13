// Searches for k such that k.100#+1 has no low factors
//
// We sieve to 2^32 on the ARM only, then after that the calculation
// of 100#^-1 mod p is offloaded to the Epiphany, with the ARM only generating
// primes and applying the result to the sieve.

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <e-hal.h>

// Epiphany data exchange structure
// Each core reads from its 16th of the prime array
// and write to its 16th of the result array.
// Todo: Keep epiphany active all the time
//       by having second set of primes ready in shared mem
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

// 26000 low primes is to just over 300000, so
// enough to sieve to 90 billion.
#define NUM_LOW_PRIMES 26000
unsigned int low_primes[NUM_LOW_PRIMES];
long long lastp = 0xffffffff;

// Main sieve size - the range of k's to consider.
//#define SIEVE_SIZE (256u*1024u*1024u*8u)
#define SIEVE_SIZE (256u*2048u*8u)
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

// The prime generator sieve for primes over 2^32 to send to the epiphany
// A worker thread generates primes by sieving, the main thread handles 
// communication with epiphany and filling results into the main sieve.
// Keep this buffer quite small to fit in the ARM L2 cache
#define GEN_PRIME_SIEVE_SIZE 64*1024*4
static unsigned int* genprimsieve;
static unsigned int* genprimsievenext;
static unsigned int genprimsievenextready;
static unsigned int genprimsievewait;
static pthread_mutex_t genprimsieve_mutex;
static pthread_cond_t genprimsieve_cv;
static pthread_cond_t genprimsievewait_cv;

static unsigned int pattern[15015];
static unsigned genprimsievej;
static unsigned int genprimoffsets[NUM_LOW_PRIMES];
static long long base;

void initpattern()
{
  for (int i = 0; i < 5; ++i)
  {
    unsigned offset = 0;
    while (offset < sizeof(pattern) * 8)
    {
      pattern[offset >> 5] |= 1<<(offset&0x1f); 
      offset += low_primes[i];
    }
  }
}

void fillsieve();
void* sievethread(void* unused);

void initgenprimes()
{
  base = lastp+2;
  int i, j;
  for (i = 5; i < NUM_LOW_PRIMES; ++i)
  {
    unsigned offset = base % low_primes[i];
    if (offset) offset = low_primes[i] - offset;
    if (offset & 1) offset += low_primes[i];
    offset >>= 1;
    genprimoffsets[i] = offset;
  }

  genprimsieve = malloc(GEN_PRIME_SIEVE_SIZE);
  genprimsievenext = malloc(GEN_PRIME_SIEVE_SIZE);
  
  fillsieve(genprimsieve, base);
  genprimsievej = 0;

  pthread_t genthread;
  genprimsievenextready = 0;
  genprimsievewait = 0;
  pthread_mutex_init(&genprimsieve_mutex, NULL);
  pthread_cond_init (&genprimsieve_cv, NULL);
  pthread_create(&genthread, NULL, sievethread, NULL);
} 

void fillsieve(unsigned int* mysieve, long long mybase)
{
  int i;
  unsigned max_prime=sqrt(mybase + GEN_PRIME_SIEVE_SIZE*16)+1;
  
  unsigned offset = mybase % 15015;
  if (offset != 0) {
    if (offset & 1) offset += 15015;
    offset >>= 1;
    offset += 15015*(offset * 9);
    offset >>= 5;
    offset %= 15015;
  }
  i = 15015 - offset;
  memcpy(mysieve, &pattern[offset], sizeof(int) * (i));
  for (; i + 15015 < GEN_PRIME_SIEVE_SIZE / sizeof(int); i+=15015)
  {
    memcpy(&mysieve[i], pattern, 15015 * sizeof(int));
  }
  memcpy(&mysieve[i], pattern, sizeof(int) * (GEN_PRIME_SIEVE_SIZE / sizeof(int) - i));
  for (i = 5; low_primes[i] < max_prime; ++i)
  {
    unsigned offset = genprimoffsets[i];
    while (offset < GEN_PRIME_SIEVE_SIZE * 8)
    {
      mysieve[offset >> 5] |= 1<<(offset&0x1f); 
      offset += low_primes[i];
    }
    genprimoffsets[i] = offset - (GEN_PRIME_SIEVE_SIZE * 8);
  }
  for (; i < NUM_LOW_PRIMES; ++i)
  {
    unsigned offset = genprimoffsets[i];
    unsigned mult = (GEN_PRIME_SIEVE_SIZE * 8) / low_primes[i];
    offset += mult * low_primes[i];
    if (offset < GEN_PRIME_SIEVE_SIZE * 8)
      offset += low_primes[i];
    genprimoffsets[i] = offset - (GEN_PRIME_SIEVE_SIZE * 8);
  }
}

void* sievethread(void* unused)
{
  while (1)
  {
    long long nextbase = base + GEN_PRIME_SIEVE_SIZE * 16;
    fillsieve(genprimsievenext, nextbase);
    
    pthread_mutex_lock(&genprimsieve_mutex);
    genprimsievenextready = 1;
    if (genprimsievewait == 1)
    {
      genprimsievewait = 0;
      pthread_cond_signal(&genprimsievewait_cv);
    }
    pthread_cond_wait(&genprimsieve_cv, &genprimsieve_mutex);
    pthread_mutex_unlock(&genprimsieve_mutex);
  }
}

// Fill buf with the next primes
void genprimes(data_t* bufptr)
{
  int i, j;

  for (i = 0; i < ENTRIES; ++genprimsievej)
  {
    if (genprimsievej == GEN_PRIME_SIEVE_SIZE * 8)
    {
      pthread_mutex_lock(&genprimsieve_mutex);
      if (!genprimsievenextready)
      {
        //fprintf(stderr, "Wait for sieve thread\n");
        genprimsievewait = 1;
        pthread_cond_wait(&genprimsievewait_cv, &genprimsieve_mutex);
      }

      base += GEN_PRIME_SIEVE_SIZE * 16;
      genprimsievej = 0;
      unsigned int* tmp = genprimsieve;
      genprimsieve = genprimsievenext;
      genprimsievenext = tmp;
      genprimsievenextready = 0;
      pthread_cond_signal(&genprimsieve_cv);
      pthread_mutex_unlock(&genprimsieve_mutex);
    }
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

// Eliminate low prime from main sieve
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

// Initialize main sieve, sieving to 2^32.
// This could be multithreaded.
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
    unsigned offset = base % 15015;
    if (offset != 0) {
      if (offset & 1) offset += 15015;
      offset >>= 1;
      offset += 15015*(offset * 9);
      offset >>= 5;
      offset %= 15015;
    }
    i = 15015 - offset;
    memcpy(sieve, &pattern[offset], sizeof(int) * (i));
    for (; i + 15015 < INIT_SIEVE_MB*1024*1024 / sizeof(int); i+=15015)
    {
      memcpy(&sieve[i], pattern, 15015 * sizeof(int));
    }
    memcpy(&sieve[i], pattern, sizeof(int) * ((INIT_SIEVE_MB*1024*1024 / sizeof(int)) - i));

    unsigned max_prime=sqrt(base + INIT_SIEVE_MB*1024*1024*16.0)+1;
    fprintf(stderr, "...\r");
    for (i = 5; low_primes[i] < max_prime; ++i)
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
//        if ((((base + (i << 1)) % 5) == 0) ||
//            (((base + (i << 1)) % 17) == 0))
//          fprintf(stderr, "Sieve broken! %u\n", base + (i << 1));
      }
      if ((i&0xffffff) == 0) 
      {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        fprintf(stderr, "%u %ld.%03ld\n", base+(i<<1), ts.tv_sec, ts.tv_nsec/1000000);
      }
    }
    if (base+(i<<1) < base)
      break;
    base += i << 1;
    fprintf(stderr, "%u\n", base);
  }
  free(sieve);
}

// Apply results from epiphany to main sieve.
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
        //if ((mainsieve[r>>5] & (1<<(r&0x1f))) == 0) printf("%lld %lld\n", r, bufptr->p[i]);
        mainsieve[r >> 5] |= 1<<(r&0x1f);
      }
    }
  }
}

int main(int argc, char*argv[])
{
  long long time = 0;
  int max_loops = 100;
  unsigned row, col, coreid, first;
  e_platform_t platform;
  e_epiphany_t dev;
  e_mem_t emem;
  data_t *curbuf = &buf, *nextbuf = &buf2;
  struct timespec sleeptime;
  sleeptime.tv_sec = 0;
  sleeptime.tv_nsec = 100000;
  
  // Generate low primes for sieving
  fprintf(stderr, "Gen Low Primes\n");
  genlowprimes();
  initpattern();

  // Allocate and init main sieve
  fprintf(stderr, "Init main sieve\n");
  mainsieve = malloc(SIEVE_SIZE >> 3);
  memset(mainsieve, 0, SIEVE_SIZE >> 3);
  initsieve();

  // Start up epiphany
  fprintf(stderr, "Init epiphany\n");
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
  fprintf(stderr, "Init gen primes\n");
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

    // Apply last results while epiphany running (if not first time)
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

    // Generate next batch of primes to send
    genprimes(nextbuf);

    // Wait for epiphany cores to complete.
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

    // Read results back
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

  // Process last set of results.
  applysieve(curbuf);

#if 0
  // Print k's with no low factor.
  for (unsigned i = 0; i < SIEVE_SIZE; ++i)
  {
    if ((mainsieve[i>>5] & (1<<(i&0x1f))) == 0)
      printf("%u\n", i);
  }
#endif

  fprintf(stderr, "Time in main loop: %lld.%03lld\n", time / 1000000000, (time / 1000000) % 1000);
  time /= max_loops;
  fprintf(stderr, "Time per loop: %lld.%04lld\n", time / 1000000000, (time / 100000) % 10000);

  e_close(&dev);
  e_free(&emem);
  e_finalize();

  return 0;
}
