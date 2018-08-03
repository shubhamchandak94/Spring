/*

BCM - A BWT-based file compressor

Copyright (C) 2008-2018 Ilya Muravyov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef _MSC_VER
#  define _FILE_OFFSET_BITS 64

#  define _fseeki64 fseeko
#  define _ftelli64 ftello
#  define _stati64 stat

#  ifdef HAVE_GETC_UNLOCKED
#    undef getc
#    define getc getc_unlocked
#  endif
#  ifdef HAVE_PUTC_UNLOCKED
#    undef putc
#    define putc putc_unlocked
#  endif
#endif

#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_DISABLE_PERFCRIT_LOCKS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "divsufsort.h" // libdivsufsort-lite

namespace spring {
namespace bcm {

typedef unsigned char byte;
typedef unsigned short word;
typedef unsigned int uint;
typedef unsigned long long ulonglong;

class bcm_class {

const char magic[5]="BCM!";

FILE* fin;
FILE* fout;

struct Encoder
{
  uint low;
  uint high;
  uint code;

  Encoder()
  {
    low=0;
    high=uint(-1);
    code=0;
  }

  void EncodeBit0(uint p, FILE *fout)
  {
#ifdef _WIN64
    low+=((ulonglong(high-low)*p)>>18)+1;
#else
    low+=((ulonglong(high-low)*(p<<(32-18)))>>32)+1;
#endif
    while ((low^high)<(1<<24))
    {
      putc(low>>24, fout);
      low<<=8;
      high=(high<<8)+255;
    }
  }

  void EncodeBit1(uint p, FILE *fout)
  {
#ifdef _WIN64
    high=low+((ulonglong(high-low)*p)>>18);
#else
    high=low+((ulonglong(high-low)*(p<<(32-18)))>>32);
#endif
    while ((low^high)<(1<<24))
    {
      putc(low>>24, fout);
      low<<=8;
      high=(high<<8)+255;
    }
  }

  void Flush(FILE *fout)
  {
    for (int i=0; i<4; ++i)
    {
      putc(low>>24, fout);
      low<<=8;
    }
  }

  void Init(FILE *fin)
  {
    for (int i=0; i<4; ++i)
      code=(code<<8)+getc(fin);
  }

  int DecodeBit(uint p, FILE *fin)
  {
#ifdef _WIN64
    const uint mid=low+((ulonglong(high-low)*p)>>18);
#else
    const uint mid=low+((ulonglong(high-low)*(p<<(32-18)))>>32);
#endif
    const int bit=(code<=mid);
    if (bit)
      high=mid;
    else
      low=mid+1;

    while ((low^high)<(1<<24))
    {
      low<<=8;
      high=(high<<8)+255;
      code=(code<<8)+getc(fin);
    }

    return bit;
  }
};

template<int RATE>
struct Counter
{
  word p;

  Counter()
  {
    p=1<<15;
  }

  void UpdateBit0()
  {
    p-=p>>RATE;
  }

  void UpdateBit1()
  {
    p+=(p^65535)>>RATE;
  }
};

struct CM: Encoder
{
  Counter<2> counter0[256];
  Counter<4> counter1[256][256];
  Counter<6> counter2[2][256][17];
  int c1;
  int c2;
  int run;

  CM()
  {
    c1=0;
    c2=0;
    run=0;

    for (int i=0; i<2; ++i)
    {
      for (int j=0; j<256; ++j)
      {
        for (int k=0; k<17; ++k)
          counter2[i][j][k].p=(k<<12)-(k==16);
      }
    }
  }

  void Encode32(uint n, FILE *fout)
  {
    for (int i=0; i<32; ++i)
    {
      if (n&(1<<31))
        Encoder::EncodeBit1(1<<17, fout);
      else
        Encoder::EncodeBit0(1<<17, fout);
      n+=n;
    }
  }

  uint Decode32(FILE *fin)
  {
    uint n=0;
    for (int i=0; i<32; ++i)
      n+=n+Encoder::DecodeBit(1<<17, fin);

    return n;
  }

  void Encode(int c, FILE *fout)
  {
    if (c1==c2)
      ++run;
    else
      run=0;
    const int f=(run>2);

    int ctx=1;
    while (ctx<256)
    {
      const int p0=counter0[ctx].p;
      const int p1=counter1[c1][ctx].p;
      const int p2=counter1[c2][ctx].p;
      const int p=((p0+p1)*7+p2+p2)>>4;

      const int j=p>>12;
      const int x1=counter2[f][ctx][j].p;
      const int x2=counter2[f][ctx][j+1].p;
      const int ssep=x1+(((x2-x1)*(p&4095))>>12);

      const int bit=c&128;
      c+=c;

      if (bit)
      {
        Encoder::EncodeBit1(ssep*3+p, fout);
        counter0[ctx].UpdateBit1();
        counter1[c1][ctx].UpdateBit1();
        counter2[f][ctx][j].UpdateBit1();
        counter2[f][ctx][j+1].UpdateBit1();
        ctx+=ctx+1;
      }
      else
      {
        Encoder::EncodeBit0(ssep*3+p,fout);
        counter0[ctx].UpdateBit0();
        counter1[c1][ctx].UpdateBit0();
        counter2[f][ctx][j].UpdateBit0();
        counter2[f][ctx][j+1].UpdateBit0();
        ctx+=ctx;
      }
    }

    c2=c1;
    c1=ctx&255;
  }

  int Decode(FILE *fin)
  {
    if (c1==c2)
      ++run;
    else
      run=0;
    const int f=(run>2);

    int ctx=1;
    while (ctx<256)
    {
      const int p0=counter0[ctx].p;
      const int p1=counter1[c1][ctx].p;
      const int p2=counter1[c2][ctx].p;
      const int p=((p0+p1)*7+p2+p2)>>4;

      const int j=p>>12;
      const int x1=counter2[f][ctx][j].p;
      const int x2=counter2[f][ctx][j+1].p;
      const int ssep=x1+(((x2-x1)*(p&4095))>>12);

      const int bit=Encoder::DecodeBit(ssep*3+p, fin);

      if (bit)
      {
        counter0[ctx].UpdateBit1();
        counter1[c1][ctx].UpdateBit1();
        counter2[f][ctx][j].UpdateBit1();
        counter2[f][ctx][j+1].UpdateBit1();
        ctx+=ctx+1;
      }
      else
      {
        counter0[ctx].UpdateBit0();
        counter1[c1][ctx].UpdateBit0();
        counter2[f][ctx][j].UpdateBit0();
        counter2[f][ctx][j+1].UpdateBit0();
        ctx+=ctx;
      }
    }

    c2=c1;
    return c1=ctx&255;
  }
} cm;

template<typename T>
inline T* mem_alloc(size_t n)
{
  T* p=(T*)malloc(n*sizeof(T));
  if (!p)
  {
    fprintf(stderr, "Out of memory!\n");
    exit(1);
  }
  return p;
}

#define mem_free(p) free(p)

void compress(int bsize)
{
  _fseeki64(fin, 0, SEEK_END);
  const long long flen=_ftelli64(fin);
  _fseeki64(fin, 0, SEEK_SET);
  if (flen>0 && bsize>flen)
    bsize=int(flen);

  byte* buf=mem_alloc<byte>(bsize);
  int* ptr=mem_alloc<int>(bsize);

  putc(magic[0], fout);
  putc(magic[1], fout);
  putc(magic[2], fout);
  putc(magic[3], fout);

  int n;
  while ((n=fread(buf, 1, bsize, fin))>0)
  {
    const int idx=divbwt(buf, buf, ptr, n);
    if (idx<1)
    {
      perror("Divbwt failed");
      exit(1);
    }

    cm.Encode32(n,fout);
	
    cm.Encode32(idx,fout);

    for (int i=0; i<n; ++i)
      cm.Encode(buf[i],fout);
  }

  cm.Encode32(0,fout); // EOF

  cm.Flush(fout);

  mem_free(buf);
  mem_free(ptr);
}

void decompress()
{
  if (getc(fin)!=magic[0]
      ||getc(fin)!=magic[1]
      ||getc(fin)!=magic[2]
      ||getc(fin)!=magic[3])
  {
    fprintf(stderr, "Not in BCM format!\n");
    exit(1);
  }

  cm.Init(fin);

  int bsize=0;
  byte* buf=NULL;
  uint* ptr=NULL;

  int n;
  while ((n=cm.Decode32(fin))>0)
  {
    if (!bsize)
    {
      if ((bsize=n)>=(1<<24)) // 5n
        buf=mem_alloc<byte>(bsize);

      ptr=mem_alloc<uint>(bsize);
    }

    const int idx=cm.Decode32(fin);
    if (n>bsize || idx<1 || idx>n)
    {
      fprintf(stderr, "Corrupt input!\n");
      exit(1);
    }

    // Inverse BW-transform
    if (n>=(1<<24)) // 5n
    {
      int t[257]={0};
      for (int i=0; i<n; ++i)
        ++t[(buf[i]=cm.Decode(fin))+1];
      for (int i=1; i<256; ++i)
        t[i]+=t[i-1];
      for (int i=0; i<n; ++i)
        ptr[t[buf[i]]++]=i+(i>=idx);
      for (int p=idx; p;)
      {
        p=ptr[p-1];
        const int c=buf[p-(p>=idx)];
        putc(c, fout);
      }
    }
    else // 4n
    {
      int t[257]={0};
      for (int i=0; i<n; ++i)
        ++t[(ptr[i]=cm.Decode(fin))+1];
      for (int i=1; i<256; ++i)
        t[i]+=t[i-1];
      for (int i=0; i<n; ++i)
        ptr[t[ptr[i]&255]++]|=(i+(i>=idx))<<8;
      for (int p=idx; p;)
      {
        p=ptr[p-1]>>8;
        const int c=ptr[p-(p>=idx)]&255;
        putc(c, fout);
      }
    }

  }
  mem_free(buf);
  mem_free(ptr);
}

public:
int bcm_main(const char *infile, const char *outfile, bool do_decomp = false, int bsize = (1<<27)-1)
{
//  const clock_t start=clock();

      if (bsize<1)
      {
        fprintf(stderr, "Block size is out of range!\n");
        exit(1);
      }
/*  if (argc<2)
  {
    fprintf(stderr,
        "BCM - A BWT-based file compressor, v1.30\n"
        "Copyright (C) 2008-2018 Ilya Muravyov\n"
        "\n"
        "Usage: %s [options] infile [outfile]\n"
        "\n"
        "Options:\n"
        "  -b#[k] Set block size to # MB or KB (default is 16 MB)\n"
        "  -d     Decompress\n"
        "  -f     Force overwrite of output file\n", argv[0]);
    exit(1);
  }
*/
  fin=fopen(infile, "rb");
  if (!fin)
  {
    perror(infile);
    exit(1);
  }

  fout=fopen(outfile, "wb");
  if (!fout)
  {
    perror(outfile);
    exit(1);
  }

  if (do_decomp)
  {
    decompress();
  }
  else
  {
    compress(bsize);
  }


  fclose(fin);
  fclose(fout);
/*
#ifndef NO_UTIME
  struct _stati64 sb;
  if (_stati64(argv[1], &sb))
  {
    perror("Stat failed");
    exit(1);
  }
  struct utimbuf ub;
  ub.actime=sb.st_atime;
  ub.modtime=sb.st_mtime;
  if (utime(ofname, &ub))
  {
    perror("Utime failed");
    exit(1);
  }
#endif
*/
  return 0;
}

}; //class bcm_class

} // namespace bcm
} // namespace spring
