#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <stdint.h>
#include <random>
#include <omp.h>

using namespace std;

#define RN    6
#define BR    32
#define BR_HALF (BR / 2)
#define PRINT_INTER 0

static const int Sbox[BR_HALF] = {
    0xc,0xa,0xd,0x3,0xe,0xb,0xf,0x7,
    0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6
};

static const int perm[BR] = {
    31,6,29,14,1,12,21,8,
    27,2,3,0,25,4,23,10,
    15,22,13,30,17,28,5,24,
    11,18,19,16,9,20,7,26
};

static const int RC0[41] = {
    0x0,0x0,0x1,0x3,0x7,0xf,0xf,0xf,0xe,0xd,
    0xa,0x5,0xa,0x5,0xb,0x6,0xc,0x9,0x3,0x6,
    0xd,0xb,0x7,0xe,0xd,0xb,0x6,0xd,0xa,0x4,
    0x9,0x2,0x4,0x9,0x3,0x7,0xe,0xc,0x8,0x1,0x2
};

static const int RC1[41] = {
    0x4,0xc,0xc,0xc,0xc,0xc,0x8,0x4,0x8,0x4,
    0x8,0x4,0xc,0x8,0x0,0x4,0xc,0x8,0x4,0xc,
    0xc,0x8,0x4,0xc,0x8,0x4,0x8,0x0,0x4,0x8,
    0x0,0x4,0xc,0xc,0x8,0x0,0x0,0x4,0x8,0x4,0xc
};

void sboxkey(int* state, const int* k, int r)
{
    for (int i = 0; i < BR_HALF; i++)
        state[i] = Sbox[state[i] & 0xF] ^ k[(r % 2) * 16 + i];
}

void permutation(int* state)
{
    int tmp[BR];
    for (int i = 0; i < BR; i++) tmp[i] = state[i];
    for (int i = 0; i < BR; i++) state[perm[i]] = tmp[i];
}

void inv_permutation(int* state)
{
    int tmp[BR];
    for (int i = 0; i < BR; i++) tmp[i] = state[i];
    for (int i = 0; i < BR; i++) state[i] = tmp[perm[i]];
}

void enc(const int* m, int* c, const int* k, int R)
{
    int state[BR], temp[BR_HALF];
    for (int i = 0; i < BR; i++) state[i] = m[i];

    for (int r = 0; r < R; r++)
    {
        for (int i = 0; i < BR_HALF; i++)
            temp[i] = state[2 * i];

        sboxkey(temp, k, r);

        for (int i = 0; i < BR_HALF; i++)
            state[2 * i + 1] ^= temp[i];

        state[1] ^= RC0[r];
        state[3] ^= RC1[r];

        permutation(state);
    }

    for (int i = 0; i < BR; i++) c[i] = state[i];
}

void dec(int* m, const int* c, const int* k, int R)
{
    int state[BR], temp[BR_HALF];
    for (int i = 0; i < BR; i++) state[i] = c[i];

    for (int r = R - 1; r >= 0; r--)
    {
        inv_permutation(state);

        state[1] ^= RC0[r];
        state[3] ^= RC1[r];

        for (int i = 0; i < BR_HALF; i++)
            temp[i] = state[2 * i];

        sboxkey(temp, k, r);

        for (int i = 0; i < BR_HALF; i++)
            state[2 * i + 1] ^= temp[i];
    }

    for (int i = 0; i < BR; i++) m[i] = state[i];
}

void convert_hexstr_to_statearray(const char* s, int* dx)
{
    for (int i = 0; i < 32; i++)
        dx[i] = (int)(strtol(string(1, s[i]).c_str(), nullptr, 16) & 0xF);
}

int main()
{
    char pt_delta_str[] = "00000000000000000000000000000000";
    char ct_delta_str[] = "00000000000000000000000000000000";

    pt_delta_str[16] = 'a';
    pt_delta_str[19] = 'a';
    ct_delta_str[23] = 'a';
    //ct_delta_str[15] = 'a';

    int pt_delta[32], ct_delta[32];
    convert_hexstr_to_statearray(pt_delta_str, pt_delta);
    convert_hexstr_to_statearray(ct_delta_str, ct_delta);

    const int R = 10;
    const long long nums = 1LL << 34;

    printf("Boomerang experiment: %lld queries, %d rounds\n", nums, R);
    printf("%s\n", pt_delta_str);
    printf("%s\n", ct_delta_str);

    int key[32];
    std::mt19937 key_rng(time(NULL));
    std::uniform_int_distribution<int> dist4(0, 15);
    for (int i = 0; i < 32; i++) key[i] = dist4(key_rng);

    long long successful_cnt = 0;

#pragma omp parallel reduction(+:successful_cnt)
    {
        std::mt19937 rng((unsigned)time(NULL) ^ omp_get_thread_num());
        std::uniform_int_distribution<int> dist(0, 15);

#pragma omp for schedule(static)
        for (long long ctr = 0; ctr < nums; ctr++)
        {
            int pt1[32], pt2[32], pt3[32], pt4[32];
            int ct1[32], ct2[32], ct3[32], ct4[32];

            for (int i = 0; i < 32; i++)
            {
                pt1[i] = dist(rng);
                pt2[i] = pt1[i] ^ pt_delta[i];
            }

            enc(pt1, ct1, key, R);
            enc(pt2, ct2, key, R);

            for (int i = 0; i < 32; i++)
            {
                ct3[i] = ct1[i] ^ ct_delta[i];
                ct4[i] = ct2[i] ^ ct_delta[i];
            }

            dec(pt3, ct3, key, R);
            dec(pt4, ct4, key, R);

            bool ok = true;
            for (int i = 0; i < 32; i++)
            {
                if ((pt3[i] ^ pt4[i]) != pt_delta[i])
                {
                    ok = false;
                    break;
                }
            }

            successful_cnt += ok;
        }
    }
    cout << "Catch: " << successful_cnt << endl;

    if (successful_cnt == 0)
        printf("overall: probability < 2^%.3f\n", -log2((double)nums));
    else
        printf("overall: 2^%.3f\n",
            log2((double)successful_cnt / (double)nums));

    return 0;
}
