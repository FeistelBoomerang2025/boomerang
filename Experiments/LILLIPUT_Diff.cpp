#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <omp.h>
#include <cmath>

using namespace std;

// 16 branches, 4 bit each branch (64-bit block)
#define BR 16

static const uint8_t SBOX[16] = { 0x4, 0x8, 0x7, 0x1, 0x9, 0x3, 0x2, 0xE, 0x0, 0xB, 0x6, 0xF, 0xA, 0x5, 0xD, 0xC };

static const uint8_t PI[16] = { 13, 9, 14, 8, 10, 11, 12, 15, 4, 5, 3, 1, 2, 6, 0, 7 };
static const uint8_t PI_inv[16] = { 14, 11, 12, 10, 8, 9, 13, 15, 3, 1, 4, 5, 6, 0, 2, 7 };

// round function of LILLIPUT
void lilliput_transform(int* state, int* round_key) {
    int next_state[BR];
    for (int i = 0; i < BR; i++) next_state[i] = state[i];

    // NonLinearLayer 
    for (int i = 0; i < 8; i++) {
        int input_val = state[7 - i] ^ round_key[i];
        int f_out = SBOX[input_val & 0xF];
        next_state[i + 8] ^= f_out;
    }

    // LinearLayer
    for (int i = 1; i < BR / 2; i++) {
        next_state[15] ^= state[i];
    }
    for (int i = 9; i < BR - 1; i++) {
        next_state[i] ^= state[BR / 2 - 1];
    }

    for (int i = 0; i < BR; i++) state[i] = next_state[i];
}

void lilliput_permutation(int* state) {
    int temp[BR];
    for (int i = 0; i < BR; i++) temp[PI[i]] = state[i];
    for (int i = 0; i < BR; i++) state[i] = temp[i];
}

void lilliput_permutation_inv(int* state) {
    int temp[BR];
    for (int i = 0; i < BR; i++) temp[PI_inv[i]] = state[i];
    for (int i = 0; i < BR; i++) state[i] = temp[i];
}

void enc(int* m, int* c, int** keys, int R) {
    int state[BR];
    for (int i = 0; i < BR; i++) state[i] = m[i] & 0xF;

    for (int r = 0; r < R; r++) {
        lilliput_transform(state, keys[r]);

        if (r < R - 1) lilliput_permutation(state);
    }

    for (int i = 0; i < BR; i++) c[i] = state[i];
}

void dec(int* c, int* m, int** keys, int R) {
    int state[BR];
    for (int i = 0; i < BR; i++) state[i] = c[i] & 0xF;

    for (int r = R - 1; r >= 0; r--) {
        lilliput_transform(state, keys[r]);

        if (r != 0) lilliput_permutation_inv(state);
    }

    for (int i = 0; i < BR; i++) m[i] = state[i];
}

void convert_hexstr_to_statearray(char* s, int dx[BR]) {
    for (int i = 0; i < BR; i++) {
        char h[2] = { s[i], '\0' };
        dx[BR - 1 - i] = strtol(h, nullptr, 16) & 0xF;
    }
}

void print(int* state) {
    for (int i = BR - 1; i >= 0; i--) {
        cout << hex << state[i];

        if (i == BR / 2) cout << ' ';
    }
    cout << endl;
}

int main() {

    char pt_delta_str[] = "0060000000000606";

    char ct_delta_str[] = "6066060000600000"; 

    int num_rounds = 4; 

    int pt_delta[BR], ct_delta[BR];
    convert_hexstr_to_statearray(pt_delta_str, pt_delta);
    convert_hexstr_to_statearray(ct_delta_str, ct_delta);

    uint64_t nums = 1ULL << 30; 

    print(pt_delta);
    print(ct_delta);
    printf("LILLIPUT Differential experiment: 2^%.1f samples, rounds = %d\n",
        log2((double)nums), num_rounds);

    uint64_t successful_cnt = 0;

#pragma omp parallel reduction(+:successful_cnt)
    {
        std::mt19937 rng((unsigned)time(NULL) ^ omp_get_thread_num());
        std::uniform_int_distribution<int> dist(0, 15);

        // assumption of random key schedule
        int* r_keys[100];
        for (int i = 0; i < num_rounds; i++) {
            r_keys[i] = new int[8];
            for (int j = 0; j < 8; j++) r_keys[i][j] = dist(rng);
        }


#pragma omp for schedule(static)
        for (uint64_t ctr = 0; ctr < nums; ctr++) {
            int pt1[BR], pt2[BR], ct1[BR], ct2[BR];

            for (int i = 0; i < BR; i++) {
                ct1[i] = dist(rng);
                ct2[i] = ct1[i] ^ ct_delta[i];
            }

            dec(ct1, pt1, r_keys, num_rounds);
            dec(ct2, pt2, r_keys, num_rounds);

            bool ok = true;
            for (int i = 0; i < BR; i++) {
                if ((pt1[i] ^ pt2[i]) != pt_delta[i]) {
                    ok = false;
                    break;
                }
            }


            successful_cnt += ok;
        }


        for (int i = 0; i < num_rounds; i++) delete[] r_keys[i];
    }

        if (successful_cnt > 0) {
            double prob = (double)successful_cnt / (double)nums;
            printf("Overall Probability = 2^%.3f\n", log2(prob));
        }
        else {
            printf("No hits found.\n");
        }
        return 0;
}
