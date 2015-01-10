#include <stdio.h>
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
  for (i = 0; i < sizeof(low_primes) / sizeof(int); ++i)
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

int main(int argc, char*argv[])
{
  genlowprimes();
  genprimes();

  unsigned row, col, coreid;
  e_platform_t platform;
  e_epiphany_t dev;
  e_mem_t emem;
  
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

  for (int loops = 0; loops < 32; ++loops)
  {
    e_write(&emem, 0, 0, 0x0, buf, _BufSize);
    for (int i = 0; i < 16; ++i)
    {
      e_start(&dev, i>>2, i&3);
    }
    genprimes();
    for (int i = 0; i < 16; ++i)
    {
      while (1)
      {
        unsigned status;
        e_read(&dev, i>>2, i&3, E_REG_STATUS, &status, sizeof(status));
        if ((status & 1) == 0)
          break;
        fprintf(stderr, ".");
      }
    }
    fprintf(stderr, "\n");
    e_read(&emem, 0, 0, 0x0, buf2, _BufSize);
  
    for (int i = 0; i < 16; ++i)
      for (int j = 0; j < _BufEntriesPerCore; j+=255)
        printf("%lld %lld\n", buf2[i*_BufEntriesPerCore+j].p, buf2[i*_BufEntriesPerCore+j].res);
  }

  e_close(&dev);
  e_free(&emem);
  e_finalize();

  return 0;
}
