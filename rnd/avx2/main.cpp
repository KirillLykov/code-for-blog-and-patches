/**
* Impelementation of SARU random number generator using SSE 4.2.
* The very first and naive version, it is 10% slower than serial code when using g++-9.2 and 35% faster with clang-602.0.49 
* Compare this implementaion towards the standard one (taken from Afshar's paper).
* References: 
* Y. Afshar, F. Schmid, A. Pishevar, and S. Worley. Exploiting seeding of random number generators for
* efficient domain decomposition parallelization of dissipative particle dynamics. Comput. Phys. Commun., 184(4):1119–1128, 2013.
*/
#include <iostream>
#include <vector>
#include <iomanip>
#include <bitset>
#include <cmath>
#include <cassert>
#include "../timer.hpp"
#include "immintrin.h"

#define noinline __attribute__((noinline))
#define finline  __attribute__((always_inline))

typedef float real;

#ifdef __AVX__
#define ALIGN alignas(32) 
#define BS 8
#define SIMD_TYPE "AVX "
#endif

#ifndef NDEBUG
// for debug
#include <cstdio>

template<typename T>
void ps(const char* s, T i)
{
    printf("%s", s);
    for (int c = 8 * sizeof(T) - 1; c >= 0; c--)
    {
        char bit = ((i >> c) & 1) ? '1' : '0';
        printf("%c", bit);
        if (c % 8 == 0) printf(" ");
    }
    printf("\n");
}

void pi(const char* s, __m256i& x, int n = 0)
{
    ALIGN int32_t buf[BS];
    _mm256_store_si256((__m256i*)&buf[0], x);
    ps(s, buf[n]);
}

void pf(const char* s, const __m256& x, int n = 0)
{
    ALIGN float buf[BS];
    _mm256_store_ps(&buf[0], x);
    ps(s, *(int32_t*)&buf[n]);
}
#endif

// 1. serial
real saruSerial(unsigned int seed1, unsigned int seed2, unsigned int seed3)
{
    seed3 ^= (seed1<<7)^(seed2>>6);
    seed2 += (seed1>>4)^(seed3>>15);
    seed1 ^= (seed2<<9)+(seed3<<8);
    
    seed3 ^= 0xA5366B4D*((seed2>>11) ^ (seed1<<1));
    seed2 += 0x72BE1579*((seed1<<4)  ^ (seed3>>16));

    seed1 ^= 0X3F38A6ED*((seed3>>5)  ^ (((signed int)seed2)>>22));
    seed2 += seed1*seed3;
    seed1 += seed3 ^ (seed2>>2);
   
    seed2 ^= ((signed int)seed2)>>17;

    int state  = 0x79dedea3*(seed1^(((signed int)seed1)>>14));
    int wstate = (state + seed2) ^ (((signed int)state)>>8);
    state  = state + (wstate*(wstate^0xdddf97f5));
    wstate = 0xABCB96F7 + (wstate>>1);
    state  = 0x4beb5d59*state + 0x2600e1f7; // LCG
    
    wstate = wstate + 0x8009d14b + ((((signed int)wstate)>>31)&0xda879add); // OWS
    unsigned int v = (state ^ (state>>26))+wstate;
    unsigned int r = (v^(v>>20))*0x6957f5a7;
    real res = r / (4294967295.0f);
    return res;
}

real randSerial(size_t i, size_t j, size_t idtimestep)
{
  const real rnd = saruSerial(std::min(i, j), std::max(i, j), idtimestep);
  return 3.464101615f * rnd - 1.732050807f;
}

#ifdef __AVX__

__m256 _mm256_cvtepu32_ps(const __m256i& v)
{
    __m256i v2 = _mm256_srli_epi32(v, 1);     // v2 = v / 2
    __m256i v1 = _mm256_sub_epi32(v, v2);     // v1 = v - (v / 2)
    __m256 v2f = _mm256_cvtepi32_ps(v2);
    __m256 v1f = _mm256_cvtepi32_ps(v1);
    return _mm256_add_ps(v2f, v1f);
}

// 2. SSE4.2
__m256 saruSSE(__m256i seed1, __m256i seed2, __m256i seed3)
{
    //seed3 ^= (seed1<<7)^(seed2>>6);
    __m256i t1 = _mm256_slli_epi32(seed1, 7);
    __m256i t2 = _mm256_srli_epi32(seed2, 6);
    __m256i t3 = _mm256_xor_si256(t1, t2);
    seed3 = _mm256_xor_si256(seed3, t3);

    //seed2 += (seed1>>4)^(seed3>>15);
    t1 = _mm256_srli_epi32(seed1, 4);
    t2 = _mm256_srli_epi32(seed2, 15);
    t3 = _mm256_xor_si256(t1, t2);
    seed2 = _mm256_add_epi32(seed2, t3);

    //seed1 ^= (seed2<<9)+(seed3<<8);
    t1 = _mm256_slli_epi32(seed2, 9);
    t2 = _mm256_slli_epi32(seed3, 8);
    t3 = _mm256_add_epi32(t1, t2);
    seed1 = _mm256_xor_si256(seed1, t3);    

    //seed3 ^= 0xA5366B4D*((seed2>>11) ^ (seed1<<1));
    t1 = _mm256_srli_epi32(seed2, 11);
    t2 = _mm256_slli_epi32(seed1, 1);
    t3 = _mm256_xor_si256(t1, t2);
    t2 = _mm256_set1_epi32(0xA5366B4D);
    t3 = _mm256_mullo_epi32(t3, t2);
    seed3 = _mm256_xor_si256(seed3, t3);

    //seed2 += 0x72BE1579*((seed1<<4)  ^ (seed3>>16));
    t1 = _mm256_slli_epi32(seed1, 4);
    t2 = _mm256_srli_epi32(seed3, 16);
    t3 = _mm256_xor_si256(t1, t2);
    t2 = _mm256_set1_epi32(0x72BE1579);
    t3 = _mm256_mullo_epi32(t2, t3);
    seed2 = _mm256_add_epi32(seed2, t3);    

    //seed1 ^= 0X3F38A6ED*((seed3>>5)  ^ (((signed int)seed2)>>22));
    t1 = _mm256_srli_epi32(seed3, 5);
    t2 = _mm256_srai_epi32(seed2, 22);
    t3 = _mm256_xor_si256(t1, t2);
    t2 = _mm256_set1_epi32(0X3F38A6ED);
    t3 = _mm256_mullo_epi32(t3, t2);
    seed1 = _mm256_xor_si256(seed1, t3);    

    //seed2 += seed1*seed3;
    t1 = _mm256_mullo_epi32(seed1, seed3);
    seed2 = _mm256_add_epi32(seed2, t1);    

    //seed1 += seed3 ^ (seed2>>2);
    t1 = _mm256_srli_epi32(seed2, 2);
    t2 = _mm256_xor_si256(seed3, t1);
    seed1 = _mm256_add_epi32(seed1, t2);

    // seed2 ^= ((signed int)seed2)>>17;
    t1 =  _mm256_srai_epi32(seed2, 17);
    seed2 = _mm256_xor_si256(seed2, t1);

    //int state  = 0x79dedea3*(seed1^(((signed int)seed1)>>14));
    t1 = _mm256_srai_epi32(seed1, 14);
    t2 = _mm256_xor_si256(seed1, t1);
    t3 = _mm256_set1_epi32(0x79dedea3);
    __m256i state = _mm256_mullo_epi32(t3, t2);
    
    //int wstate = (state + seed2) ^ (((signed int)state)>>8);
    t1 = _mm256_srai_epi32(state, 8);
    t2 = _mm256_add_epi32(state, seed2);
    __m256i wstate = _mm256_xor_si256(t1, t2);

    //state  = state + (wstate*(wstate^0xdddf97f5));
    t1 = _mm256_set1_epi32(0xdddf97f5);    
    t2 = _mm256_xor_si256(wstate, t1);
    t3 = _mm256_mullo_epi32(wstate, t2);
    state = _mm256_add_epi32(state, t3);

    //wstate = 0xABCB96F7 + (wstate>>1);
    t1 = _mm256_srai_epi32(wstate, 1);
    t2 = _mm256_set1_epi32(0xABCB96F7);
    wstate = _mm256_add_epi32(t2, t1);
    //state  = 0x4beb5d59*state + 0x2600e1f7; // LCG
    t1 = _mm256_set1_epi32(0x4beb5d59);    
    t2 = _mm256_mullo_epi32(t1, state);
    t3 = _mm256_set1_epi32(0x2600e1f7);
    state = _mm256_add_epi32(t2, t3);

    //wstate = wstate + 0x8009d14b + ((((signed int)wstate)>>31)&0xda879add); // OWS
    t1 = _mm256_srai_epi32(wstate, 31);
    t2 = _mm256_set1_epi32(0xda879add);
    t3 = _mm256_and_si256(t1, t2);
    t1 = _mm256_set1_epi32(0x8009d14b);    
    t2 = _mm256_add_epi32(t1, t3);
    wstate = _mm256_add_epi32(wstate, t2);

    //unsigned int v = (state ^ (state>>26))+wstate;
    t1 = _mm256_srai_epi32(state, 26);
    t2 = _mm256_xor_si256(state, t1);   
    t3 = _mm256_add_epi32(t2, wstate);
    //unsigned int r = (v^(v>>20))*0x6957f5a7;
    t1 = _mm256_srli_epi32(t3, 20);
    t2 = _mm256_xor_si256(t3, t1);
    t3 = _mm256_set1_epi32(0x6957f5a7);
    t1 = _mm256_mullo_epi32(t2, t3);
    //real res = r / (4294967295.0);
    __m256 res = _mm256_cvtepu32_ps(t1);
    __m256 scal = _mm256_set1_ps(4294967295.0f);
    return _mm256_div_ps(res, scal);
}

// all arrays must be alignas(32)
void randSSE(uint32_t* i, uint32_t* j, uint32_t* idtimestep, float* res)
{   
    __m256i ti = _mm256_load_si256((__m256i *)i);
    __m256i tj = _mm256_load_si256((__m256i *)j);
    __m256i tidts = _mm256_load_si256((__m256i *)idtimestep);
    __m256i tmin = _mm256_min_epu32(ti, tj);
    __m256i tmax = _mm256_max_epu32(ti, tj);

    __m256 rnd = saruSSE(tmin, tmax, tidts);
    __m256 scal1 = _mm256_set1_ps(3.464101615f);
    rnd = _mm256_mul_ps(rnd, scal1); 
    scal1 = _mm256_set1_ps(1.732050807f);
    rnd = _mm256_sub_ps(rnd, scal1);
    _mm256_store_ps(res, rnd);    
}
#endif

void checkCorrectness()
{
    const unsigned int ngenerate = BS*10000;
    float data[ngenerate];
    ALIGN unsigned int seed1[ngenerate];
    ALIGN unsigned int seed2[ngenerate];
    ALIGN unsigned int seed3[ngenerate];
    for (unsigned int i = 0; i < ngenerate; ++i) {
        seed1[i] = (rand()+1)%10;
        seed2[i] = (rand()+2)%100;
        seed3[i] = (rand()+3)%1000;
    }
     
    for (int i = 0; i < ngenerate; ++i) {
        data[i] = randSerial(seed1[i], seed2[i], seed3[i]);
    }

    ALIGN float dataSSE[ngenerate];
    for (int i = 0; i < ngenerate; i += BS) {
        randSSE(&seed1[i], &seed2[i], &seed3[i], &dataSSE[i]);
    }

    // correctness
    float eps = 1e-6;
    for (unsigned int i = 0; i < ngenerate; ++i) {
        assert( fabs(dataSSE[i] - data[i]) < eps );     
    }
} 

void checkSpeed()
{
    const unsigned int ngenerate = BS*1000;
    ALIGN unsigned int seed[3*ngenerate];
    for (unsigned int i = 0; i < 3*ngenerate; ++i) {
        seed[i] = (rand()+i%3)%(100);
    }
 
    {   float tmp = 0.0f; // to prevent optimizing out
        timer t;
        t.start();
        for (int nruns = 0; nruns < 10000; ++nruns)
            for (int i = 0; i < ngenerate; i+=3) {
                tmp += randSerial(seed[i], seed[i+1], seed[i+2]);
            }
        uint64_t tim = t.stop();
        std::cout << "Ser " << tim << ". " << tmp << std::endl;  
    }
    {
        timer t;
        t.start();
        ALIGN float buf[BS];
        float tmp = 0.0f;
        for (int nruns = 0; nruns < 10000; ++nruns)
            for (int i = 0; i < ngenerate; i += 3*BS) {
                randSSE(&seed[i], &seed[i + BS], &seed[i + 2*BS], &buf[0]);
                tmp += buf[0] + buf[1] + buf[2] + buf[3] + buf[4] + buf[5] +buf[6] + buf[7]; // do it to prevent compiler remove unused 
            }
        uint64_t tim = t.stop();
        std::cout << SIMD_TYPE << tim << ". "<< tmp << std::endl;
    }
}

void checkSimple()
{
    // an example of loading constants
    //__m256 first = _mm256_set_ss(1.732050807f);
    //first = _mm256_shuffle_ps(first, first, 0x00);
    //__m256i second = _mm256_cvtsi32_si128(0X3F38A6ED);
    //second = _mm256_shuffle_epi32(second, 0x00);
    //pi("A1 ", second, 0);
    //pi("A2 ", second, 1);
    //pi("A3 ", second, 2);
    //pi("A4 ", second, 3);   
 
    ALIGN unsigned int buf1[8] = {2, 1, 3, 5, 5, 6, 7, 8};
    ALIGN unsigned int buf2[8] = {3, 2, 4, 4, 5, 6, 7, 8};
    ALIGN unsigned int buf3[8] = {4, 3, 5, 6, 5, 6, 7, 8};
    ALIGN float res[8];
    randSSE(buf1, buf2, buf3, res);

    for (int i = 0; i < 8; ++i)
        std::cout << "Orig = " << randSerial(buf1[i], buf2[i], buf3[i]) << ". Res = " << res[i] <<std::endl;
    
}

int main(int argc, char *argv[])
{
    checkCorrectness();
    checkSpeed();
    return 0;
}
