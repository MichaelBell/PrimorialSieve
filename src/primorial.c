#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <e-hal.h>

typedef struct inout_data
{
  long long p;
  long long res;
} data_t;

#define _BufEntriesPerCore 1024
#define _BufEntries (_BufEntriesPerCore*16)
#define _BufSize   (sizeof(buf))
#define _BufOffset (0x01000000)

data_t buf[_BufEntries];
data_t buf2[_BufEntries];
unsigned int low_primes[100000];
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

// Fill buf with the next primes
void genprimes()
{
  unsigned int sieve[5800] = {0};
  long long base = lastp+2;
  int i, j;
  unsigned max_prime=sqrt(base + sizeof(sieve)*16)+1;
  for (i = 0; low_primes[i] < max_prime; ++i)
  {
    long long offset = base % low_primes[i];
    if (offset) offset = low_primes[i] - offset;
    if (offset & 1) offset += low_primes[i];
    offset >>= 1;
    while (offset < sizeof(sieve) * 8)
    {
      sieve[offset >> 5] |= 1<<(offset&0x1f); 
      offset += low_primes[i];
    }
  }
  for (i = 0, j = 0; i < sizeof(buf) / sizeof(data_t); ++j)
  {
    if ((sieve[j>>5] & (1 << (j&0x1f))) == 0)
    {
      buf[i].p = base + (j << 1);
      buf[i++].res = 0;
    }
  }
  lastp = buf[i-1].p;
  fprintf(stderr, "Oversieved by %d\n", sizeof(sieve)*8 - j);
}

// return t such that at = 1 mod m
unsigned inverse(unsigned a, unsigned m)
{
  long long t = 0, newt = 1, r = m, newr = a;
  while (newr != 0)
  {
    long long q = r / newr;
    long long x = t - q * newt;
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

void sievep(unsigned p)
{
  unsigned n = (223092870ll * 2756205443ll) % p;
  n = (n * 907383479ll) % p;
  n = (n * 4132280413ll) % p;
  n = inverse(n, p);
  unsigned r = (n == 0) ? 0 : (p - n);
  if (p==1795133)
        fprintf(stderr, "%u %u\n", p, r);
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
      if ((i&0xffffff) == 0) fprintf(stderr, "%u\n", base+(i<<1));
    }
    if (base+(i<<1) < base)
      break;
    base += i << 1;
    fprintf(stderr, "%u\n", base);
  }
  free(sieve);
}

void applysieve()
{
  for (int i = 0; i < _BufEntries; ++i)
  {
    unsigned long long r = (buf2[i].res == 0) ? 0 : (buf2[i].p - buf2[i].res);
    if (r < SIEVE_SIZE)
    {
      if ((mainsieve[r>>5] & (1<<(r&0x1f))) == 0) printf("%lld %lld\n", r, buf2[i].p);
      mainsieve[r >> 5] |= 1<<(r&0x1f);
    }
  }
}

int main(int argc, char*argv[])
{
  unsigned row, col, coreid, first;
  e_platform_t platform;
  e_epiphany_t dev;
  e_mem_t emem;
  struct timespec sleeptime;
  sleeptime.tv_sec = 0;
  sleeptime.tv_nsec = 1000000;
  
  // Generate low primes for sieving
  genlowprimes();

  // Allocate and init main sieve
  mainsieve = malloc(SIEVE_SIZE >> 3);
  memset(mainsieve, 0, SIEVE_SIZE >> 3);
  initsieve();

  // Start up epiphany
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
  genprimes();
  first = 1;

  // Main loop
  for (int loops = 0; loops < 128; ++loops)
  {
    // Copy in next batch of primes for epiphany and start the cores
    e_write(&emem, 0, 0, 0x0, buf, _BufSize);
    for (int i = 0; i < 16; ++i)
    {
      e_start(&dev, i>>2, i&3);
    }
    
    // Generate next batch while they are running
    genprimes();
    
    // And apply last results (if not first time)
    if (!first)
      applysieve();
    else
      first = 0;

    // Wait for epiphany cores.
    for (int i = 0; i < 16; ++i)
    {
      while (1)
      {
        unsigned status;
        e_read(&dev, i>>2, i&3, E_REG_STATUS, &status, sizeof(status));
        if ((status & 1) == 0)
          break;
        fprintf(stderr, ".");
        nanosleep(&sleeptime, NULL);
      }
    }
    e_read(&emem, 0, 0, 0x0, buf2, _BufSize);
    fprintf(stderr, "\n%lld\n", buf2[_BufEntries-1].p);
  
#if 0
    for (int i = 0; i < 16; ++i)
      for (int j = 0; j < _BufEntriesPerCore; j+=255)
        printf("%lld %lld\n", buf2[i*_BufEntriesPerCore+j].p, buf2[i*_BufEntriesPerCore+j].res);
#endif
  }

  applysieve();

  for (unsigned i = 0; i < SIEVE_SIZE; ++i)
  {
    if ((mainsieve[i>>5] & (1<<(i&0x1f))) == 0)
      printf("%u\n", i);
  }

  e_close(&dev);
  e_free(&emem);
  e_finalize();

  return 0;
}
