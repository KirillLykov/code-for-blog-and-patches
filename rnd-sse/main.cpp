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
#include "timer.hpp"
#include "immintrin.h"

#define noinline __attribute__((noinline))
#define finline  __attribute__((always_inline))

typedef float real;

// for debug
#include <stdio.h>      /* printf */
#include <string.h>     /* strcat */
#include <stdlib.h>     /* strtol */

#ifndef NDEBUG
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

void pi(const char* s, __m128i& x, int n = 0)
{
    alignas(16) int32_t buf[4];
    _mm_store_si128((__m128i*)&buf[0], x);
    ps(s, buf[n]);
}

void pf(const char* s, const __m128& x, int n = 0)
{
    alignas(16) float buf[4];
    _mm_store_ps(&buf[0], x);
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

__m128 _mm_cvtepu32_ps(const __m128i& v)
{
    __m128i v2 = _mm_srli_epi32(v, 1);     // v2 = v / 2
    __m128i v1 = _mm_sub_epi32(v, v2);     // v1 = v - (v / 2)
    __m128 v2f = _mm_cvtepi32_ps(v2);
    __m128 v1f = _mm_cvtepi32_ps(v1);
    return _mm_add_ps(v2f, v1f);
}

// 2. SSE4.2
__m128 saruSSE(__m128i seed1, __m128i seed2, __m128i seed3)
{
    //seed3 ^= (seed1<<7)^(seed2>>6);
    __m128i t1 = _mm_slli_epi32(seed1, 7);
    __m128i t2 = _mm_srli_epi32(seed2, 6);
    __m128i t3 = _mm_xor_si128(t1, t2);
    seed3 = _mm_xor_si128(seed3, t3);

    //seed2 += (seed1>>4)^(seed3>>15);
    t1 = _mm_srli_epi32(seed1, 4);
    t2 = _mm_srli_epi32(seed2, 15);
    t3 = _mm_xor_si128(t1, t2);
    seed2 = _mm_add_epi32(seed2, t3);

    //seed1 ^= (seed2<<9)+(seed3<<8);
    t1 = _mm_slli_epi32(seed2, 9);
    t2 = _mm_slli_epi32(seed3, 8);
    t3 = _mm_add_epi32(t1, t2);
    seed1 = _mm_xor_si128(seed1, t3);    

    //seed3 ^= 0xA5366B4D*((seed2>>11) ^ (seed1<<1));
    t1 = _mm_srli_epi32(seed2, 11);
    t2 = _mm_slli_epi32(seed1, 1);
    t3 = _mm_xor_si128(t1, t2);
    t2 = _mm_set1_epi32(0xA5366B4D);
    t3 = _mm_mullo_epi32(t3, t2);
    seed3 = _mm_xor_si128(seed3, t3);

    //seed2 += 0x72BE1579*((seed1<<4)  ^ (seed3>>16));
    t1 = _mm_slli_epi32(seed1, 4);
    t2 = _mm_srli_epi32(seed3, 16);
    t3 = _mm_xor_si128(t1, t2);
    t2 = _mm_set1_epi32(0x72BE1579);
    t3 = _mm_mullo_epi32(t2, t3);
    seed2 = _mm_add_epi32(seed2, t3);    

    //seed1 ^= 0X3F38A6ED*((seed3>>5)  ^ (((signed int)seed2)>>22));
    t1 = _mm_srli_epi32(seed3, 5);
    t2 = _mm_srai_epi32(seed2, 22);
    t3 = _mm_xor_si128(t1, t2);
    t2 = _mm_set1_epi32(0X3F38A6ED);
    t3 = _mm_mullo_epi32(t3, t2);
    seed1 = _mm_xor_si128(seed1, t3);    

    //seed2 += seed1*seed3;
    t1 = _mm_mullo_epi32(seed1, seed3);
    seed2 = _mm_add_epi32(seed2, t1);    

    //seed1 += seed3 ^ (seed2>>2);
    t1 = _mm_srli_epi32(seed2, 2);
    t2 = _mm_xor_si128(seed3, t1);
    seed1 = _mm_add_epi32(seed1, t2);

    // seed2 ^= ((signed int)seed2)>>17;
    t1 =  _mm_srai_epi32(seed2, 17);
    seed2 = _mm_xor_si128(seed2, t1);

    //int state  = 0x79dedea3*(seed1^(((signed int)seed1)>>14));
    t1 = _mm_srai_epi32(seed1, 14);
    t2 = _mm_xor_si128(seed1, t1);
    t3 = _mm_set1_epi32(0x79dedea3);
    __m128i state = _mm_mullo_epi32(t3, t2);
    
    //int wstate = (state + seed2) ^ (((signed int)state)>>8);
    t1 = _mm_srai_epi32(state, 8);
    t2 = _mm_add_epi32(state, seed2);
    __m128i wstate = _mm_xor_si128(t1, t2);

    //state  = state + (wstate*(wstate^0xdddf97f5));
    t1 = _mm_set1_epi32(0xdddf97f5);    
    t2 = _mm_xor_si128(wstate, t1);
    t3 = _mm_mullo_epi32(wstate, t2);
    state = _mm_add_epi32(state, t3);

    //wstate = 0xABCB96F7 + (wstate>>1);
    t1 = _mm_srai_epi32(wstate, 1);
    t2 = _mm_set1_epi32(0xABCB96F7);
    wstate = _mm_add_epi32(t2, t1);
    //state  = 0x4beb5d59*state + 0x2600e1f7; // LCG
    t1 = _mm_set1_epi32(0x4beb5d59);    
    t2 = _mm_mullo_epi32(t1, state);
    t3 = _mm_set1_epi32(0x2600e1f7);
    state = _mm_add_epi32(t2, t3);

    //wstate = wstate + 0x8009d14b + ((((signed int)wstate)>>31)&0xda879add); // OWS
    t1 = _mm_srai_epi32(wstate, 31);
    t2 = _mm_set1_epi32(0xda879add);
    t3 = _mm_and_si128(t1, t2);
    t1 = _mm_set1_epi32(0x8009d14b);    
    t2 = _mm_add_epi32(t1, t3);
    wstate = _mm_add_epi32(wstate, t2);

    //unsigned int v = (state ^ (state>>26))+wstate;
    t1 = _mm_srai_epi32(state, 26);
    t2 = _mm_xor_si128(state, t1);   
    t3 = _mm_add_epi32(t2, wstate);
    //unsigned int r = (v^(v>>20))*0x6957f5a7;
    t1 = _mm_srli_epi32(t3, 20);
    t2 = _mm_xor_si128(t3, t1);
    t3 = _mm_set1_epi32(0x6957f5a7);
    t1 = _mm_mullo_epi32(t2, t3);
    //real res = r / (4294967295.0);
    __m128 res = _mm_cvtepu32_ps(t1);
    __m128 scal = _mm_set1_ps(4294967295.0f);
    return _mm_div_ps(res, scal);
}

// all arrays must be alignas(16)
void randSSE(uint32_t* i, uint32_t* j, uint32_t* idtimestep, float* res)
{   
    __m128i ti = _mm_load_si128((__m128i *)i);
    __m128i tj = _mm_load_si128((__m128i *)j);
    __m128i tidts = _mm_load_si128((__m128i *)idtimestep);
    __m128i tmin = _mm_min_epu32(ti, tj);
    __m128i tmax = _mm_max_epu32(ti, tj);

    __m128 rnd = saruSSE(tmin, tmax, tidts);
    __m128 scal1 = _mm_set1_ps(3.464101615f);
    rnd = _mm_mul_ps(rnd, scal1); 
    scal1 = _mm_set1_ps(1.732050807f);
    rnd = _mm_sub_ps(rnd, scal1);
    _mm_store_ps(res, rnd);    
}

void checkCorrectness()
{
    const unsigned int ngenerate = 4*10000;
    float data[ngenerate];
    alignas(16) unsigned int seed1[ngenerate];
    alignas(16) unsigned int seed2[ngenerate];
    alignas(16) unsigned int seed3[ngenerate];
    for (unsigned int i = 0; i < ngenerate; ++i) {
        seed1[i] = (rand()+1)%10;
        seed2[i] = (rand()+2)%100;
        seed3[i] = (rand()+3)%1000;
    }
     
    for (int i = 0; i < ngenerate; ++i) {
        data[i] = randSerial(seed1[i], seed2[i], seed3[i]);
    }

    alignas(16) float dataSSE[ngenerate];
    for (int i = 0; i < ngenerate; i += 4) {
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
    const unsigned int ngenerate = 4*1000;
    alignas(16) unsigned int seed[3*ngenerate];
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
        alignas(16) float buf[4];
        float tmp = 0.0f;
        for (int nruns = 0; nruns < 10000; ++nruns)
            for (int i = 0; i < ngenerate; i += 12) {
                randSSE(&seed[i], &seed[i+4], &seed[i+8], &buf[0]);
                tmp += buf[0] + buf[1] + buf[2] +buf[3];
            }
        uint64_t tim = t.stop();
        std::cout << "SSE " << tim << ". "<< tmp << std::endl;
    }
}

void checkSimple()
{

    __m128 first = _mm_set_ss(1.732050807f);
    first = _mm_shuffle_ps(first, first, 0x00);
    __m128i second = _mm_cvtsi32_si128(0X3F38A6ED);
    second = _mm_shuffle_epi32(second, 0x00);
    pi("A ", second, 0);
    pi("A ", second, 1);
    pi("A ", second, 2);
    pi("A ", second, 3);   
 
    unsigned int buf1[4] = {2, 1, 3, 4};
    unsigned int buf2[4] = {3, 2, 4, 5};
    unsigned int buf3[4] = {4, 3, 5, 6};
    float res[4];
    randSSE(buf1, buf2, buf3, res);

    std::cout << "Orig = " << randSerial(2, 1, 3) << ". Res = " << res[0] <<std::endl;
    std::cout << "Orig = " << randSerial(3, 2, 4) << ". Res = " << res[1] <<std::endl;
    std::cout << "Orig = " << randSerial(4, 3, 5) << ". Res = " << res[2] <<std::endl;
    std::cout << "Orig = " << randSerial(5, 4, 6) << ". Res = " << res[3] <<std::endl;
}

int main(int argc, char *argv[])
{
    checkSimple();
    return 0;
}
