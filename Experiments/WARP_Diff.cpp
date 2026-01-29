#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <omp.h>

using namespace std;

#define BR        32
#define BR_HALF   (BR / 2)
#define PRINT_INTER 0

static const int Sbox[BR_HALF] = {
    0xc, 0xa, 0xd, 0x3, 0xe, 0xb, 0xf, 0x7,
    0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6
};

static const int perm[BR] = {
    31, 6, 29, 14, 1, 12, 21, 8,
    27, 2, 3, 0, 25, 4, 23, 10,
    15, 22, 13, 30, 17, 28, 5, 24,
    11, 18, 19, 16, 9, 20, 7, 26
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

void sboxkey(int* state, int* k, int r)
{
    for (int i = 0; i < BR_HALF; i++)
        state[i] = Sbox[state[i] & 0xF] ^ k[(r & 1) * 16 + i];
}

void permutation(int* state)
{
    int tmp[BR];
    for (int i = 0; i < BR; i++) tmp[i] = state[i];
    for (int i = 0; i < BR; i++) state[perm[i]] = tmp[i];
}

void enc(int* m, int* c, int* k, int R)
{
    int state[BR], temp[BR_HALF];
    for (int i = 0; i < BR; i++) state[i] = m[i];

    for (int r = 0; r < R; r++) {
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

void convert_hexstr_to_statearray(const char* s, int dx[32])
{
    for (int i = 0; i < 32; i++) {
        char h[2] = { s[i], '\0' };
        dx[i] = strtol(h, nullptr, 16) & 0xF;
    }
}

int main()
{
    const char pt_delta_str[] = "aa00aaaa0a000000000aaa00000a0a00";
    const char ct_delta_str[] = "0a000000000a00000000000000000000";

    int pt_delta[32], ct_delta[32];
    convert_hexstr_to_statearray(pt_delta_str, pt_delta);
    convert_hexstr_to_statearray(ct_delta_str, ct_delta);

    int num_rounds = 7;

    uint64_t nums = 1ULL << 34;

    uint64_t successful_cnt = 0;

    printf("Differential experiment: 2^%.1f samples, rounds = %d\n",
        log2((double)nums), num_rounds);

#pragma omp parallel
    {
        std::mt19937 rng(
            (unsigned)time(NULL) ^ omp_get_thread_num()
        );
        std::uniform_int_distribution<int> dist(0, 15);

#pragma omp for reduction(+:successful_cnt)
        for (uint64_t ctr = 0; ctr < nums; ctr++) {

            int key[32], pt1[32], pt2[32], ct1[32], ct2[32];

            for (int i = 0; i < 32; i++)
                key[i] = dist(rng);

            for (int i = 0; i < 32; i++) {
                pt1[i] = dist(rng);
                pt2[i] = pt1[i] ^ pt_delta[i];
            }

            enc(pt1, ct1, key, num_rounds);
            enc(pt2, ct2, key, num_rounds);

            bool ok = true;
            for (int i = 0; i < 32; i++) {
                if ((ct1[i] ^ ct2[i]) != ct_delta[i]) {
                    ok = false;
                    break;
                }
            }

            successful_cnt += ok;
        }
    }

    cout << "successful_cnt = " << successful_cnt << endl;

    double prob = (double)successful_cnt / (double)nums;
    printf("overall probability = 2^%.3f\n", log2(prob));

    return 0;
}
