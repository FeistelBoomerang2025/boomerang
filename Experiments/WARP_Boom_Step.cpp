#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/random.h>
#include <random>
#include <mutex>
#include <atomic>
using namespace std;
mutex stdout_mutex;

#define RN    6    /*rotation number*/
#define BR    32  /*brunch number*/
#define BR_HALF    (BR / 2)  /*half of the branch number*/
#define PRINT_INTER 0


static const int Sbox[BR_HALF] = { 0xc, 0xa, 0xd, 0x3, 0xe, 0xb, 0xf, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6 };
static const int perm[BR] = { 31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10, 15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26, };
static const int RC0[41] = { 0x0U, 0x0U, 0x1U, 0x3U, 0x7U, 0xfU, 0xfU, 0xfU, 0xeU, 0xdU, 0xaU, 0x5U, 0xaU, 0x5U, 0xbU, 0x6U, 0xcU, 0x9U, 0x3U, 0x6U, 0xdU, 0xbU, 0x7U, 0xeU, 0xdU, 0xbU, 0x6U, 0xdU, 0xaU, 0x4U, 0x9U, 0x2U, 0x4U, 0x9U, 0x3U, 0x7U, 0xeU, 0xcU, 0x8U, 0x1U, 0x2U };
static const int RC1[41] = { 0x4U, 0xcU, 0xcU, 0xcU, 0xcU, 0xcU, 0x8U, 0x4U, 0x8U, 0x4U, 0x8U, 0x4U, 0xcU, 0x8U, 0x0U, 0x4U, 0xcU, 0x8U, 0x4U, 0xcU, 0xcU, 0x8U, 0x4U, 0xcU, 0x8U, 0x4U, 0x8U, 0x0U, 0x4U, 0x8U, 0x0U, 0x4U, 0xcU, 0xcU, 0x8U, 0x0U, 0x0U, 0x4U, 0x8U, 0x4U, 0xcU };

void printState(int* state)
{
    printf("L: ");
    for (int x = 0; x < BR_HALF; x++)
    {
        printf("%x ", state[2 * x + 0]);
    }
    printf("R: ");
    for (int x = 0; x < BR_HALF; x++)
    {
        printf("%x ", state[2 * x + 1]);
    }
    printf("\n");
}

void sboxkey(int* state, int* k, int r)
{
    for (int i = 0; i < BR_HALF; i++)
    {
        state[i] = Sbox[state[i] & 0xF] ^ k[(r % 2) * 16 + i];
    }
}

void permutation(int* state)
{
    int temp[BR];
    for (int j = 0; j < BR; j++)
    {
        temp[j] = state[j];
    }
    for (int j = 0; j < BR; j++)
    {
        state[perm[j]] = temp[j];
    }
}

void inv_permutation(int* state)
{
    int temp[BR];
    for (int j = 0; j < BR; j++)
    {
        temp[j] = state[j];
    }
    for (int j = 0; j < BR; j++)
    {
        state[j] = temp[perm[j]];
    }
}

void enc(int* m, int* c, int* k, int R)
{
    /*intermediate value*/
    int state[BR];

    /*left half intermediate value*/
    int temp[BR_HALF];

    for (int i = 0; i < BR; i++)
    {
        state[i] = m[i];
    }

    /*round function(1 to 40 round)*/
    for (int i = 0; i < R; i++)
    {
#if PRINT_INTER
        printf("%d round\n", i + 1);
        printState(state);
#endif

        for (int j = 0; j < BR_HALF; j++)
        {
            temp[j] = state[j * 2];
        }
        /*insert key and Sbox*/
        sboxkey(temp, k, i);
        /*XOR*/
        for (int j = 0; j < BR_HALF; j++)
        {
            state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
        }
        /*add round constants*/
        state[1] = state[1] ^ RC0[i];
        state[3] = state[3] ^ RC1[i];

        /*permutation*/
        permutation(state);
    }

    /*last round function */
#if PRINT_INTER
    printf("%d round\n", R);
    printState(state);
#endif

#if PRINT_INTER
    printState(state);
#endif

    /*no permutation in the last round*/

    /*copy ciphertext*/
    for (int i = 0; i < BR; i++)
    {
        c[i] = state[i];
    }

}

void enc_Step(int (*m)[32], int (*c)[32], int* k, int R)
{
    /*intermediate value*/
    int state[BR];

    /*left half intermediate value*/
    int temp[BR_HALF];

    for (int i = 0; i < BR; i++)
    {
        state[i] = m[0][i];
    }

    /*round function(1 to 40 round)*/
    for (int i = 0; i < R; i++)
    {
        for (int j = 0; j < BR; j++) {
            m[i][j] = state[j];
        }

        for (int j = 0; j < BR_HALF; j++)
        {
            temp[j] = state[j * 2];
        }
        /*insert key and Sbox*/
        sboxkey(temp, k, i);
        /*XOR*/
        for (int j = 0; j < BR_HALF; j++)
        {
            state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
        }
        /*add round constants*/
        state[1] = state[1] ^ RC0[i];
        state[3] = state[3] ^ RC1[i];

        /*permutation*/
        permutation(state);
    }

    /*last round function */

    /*no permutation in the last round*/

    /*copy ciphertext*/
    for (int i = 0; i < BR; i++)
    {
        c[R][i] = state[i];
        m[R][i] = state[i];
    }

}

void dec(int* m, int* c, int* k, int R)
{
    /*intermediate value*/
    int state[BR];

    /*left half intermediate value*/
    int temp[BR_HALF];

    for (int i = 0; i < BR; i++)
    {
        state[i] = c[i];
    }

    /*round function(1 to 40 round)*/
    for (int i = R - 1; i >= 0; i--)
    {
#if PRINT_INTER
        printf("%d round\n", i + 1);
        printState(state);
#endif

        /*inverse of permutation*/
        inv_permutation(state);
        /*add round constants*/
        state[1] = state[1] ^ RC0[i];
        state[3] = state[3] ^ RC1[i];

        for (int j = 0; j < BR_HALF; j++)
        {
            temp[j] = state[j * 2];
        }
        /*insert key and Sbox*/
        sboxkey(temp, k, i);
        /*XOR*/
        for (int j = 0; j < BR_HALF; j++)
        {
            state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
        }
    }

    /*last round function */
#if PRINT_INTER
    printf("%d round\n", R);
    printState(state);
#endif

#if PRINT_INTER
    printState(state);
#endif

    /*no permutation in the last round*/

    /*copy plaintext*/
    for (int i = 0; i < BR; i++)
    {
        m[i] = state[i];
    }

}

void dec_Step(int (*m)[32], int (*c)[32], int* k, int R)
{
    /*intermediate value*/
    int state[BR];

    /*left half intermediate value*/
    int temp[BR_HALF];

    for (int i = 0; i < BR; i++)
    {
        state[i] = c[R][i];
    }

    /*round function(1 to 40 round)*/
    for (int i = R - 1; i >= 0; i--)
    {
        /*inverse of permutation*/
        inv_permutation(state);
        /*add round constants*/
        state[1] = state[1] ^ RC0[i];
        state[3] = state[3] ^ RC1[i];

        for (int j = 0; j < BR_HALF; j++)
        {
            temp[j] = state[j * 2];
        }
        /*insert key and Sbox*/
        sboxkey(temp, k, i);
        /*XOR*/
        for (int j = 0; j < BR_HALF; j++)
        {
            state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
        }

        for (int j = 0; j < BR; j++) {
            c[i][j] = state[j];
        }
    }

    /*last round function */

    /*no permutation in the last round*/

    /*copy plaintext*/
    for (int i = 0; i < BR; i++)
    {
        m[0][i] = state[i];
    }

}

void convert_hexstr_to_statearray(char* hex_str, int dx[32])
{
    for (int i = 0; i < 32; i++)
    {
        char hex[2];
        hex[0] = hex_str[i];
        hex[1] = '\0';
        dx[i] = (int)(strtol(hex, NULL, 16) & 0xf);
    }
}

int check_all_equal_xor(int (*pt1)[32], int (*pt2)[32], int (*pt3)[32], int (*pt4)[32], int r) {

    for (int j = 0; j < 32; j++) {
        if ((pt1[r][j] ^ pt2[r][j]) != (pt3[r][j] ^ pt4[r][j])) {
            return 0;  
        }
    }
    return 1;  
}

 // Encryption by Step Test
int main()
{
    char pt_delta_str[] = "00000000000000000000000000000000";
    char ct_delta_str[] = "00000000000000000000000000000000";
    int pt_delta[32];
    int ct_delta[32];
    int num_rounds = 11;

    int nums = 1 << 20;

    printf("analyzing 2^%.1f queries * keys of the follwoing boomerang for %d rounds\n",
        log2(nums), num_rounds);

    // int successful_cnt = 0;
    srand((unsigned)time(NULL));

    pt_delta_str[19] = 'a';
    ct_delta_str[24] = 'a';
    // ct_delta_str[19] = 'a';

    convert_hexstr_to_statearray(pt_delta_str, pt_delta);
    convert_hexstr_to_statearray(ct_delta_str, ct_delta);

    printf("%s\n", pt_delta_str);
    printf("%s\n", ct_delta_str);

    int successful_cnt = 0;
    int check_cnt = 0;


    for (int ctr = 0; ctr < nums; ++ctr)
    {
        int key[32];
        int rows = num_rounds + 1;
        int (*pt1)[32] = new int[rows][32];
        int (*pt2)[32] = new int[rows][32];
        int (*pt3)[32] = new int[rows][32];
        int (*pt4)[32] = new int[rows][32];
        int (*ct1)[32] = new int[rows][32];
        int (*ct2)[32] = new int[rows][32];
        int (*ct3)[32] = new int[rows][32];
        int (*ct4)[32] = new int[rows][32];
        for (int i = 0; i < 32; ++i) key[i] = rand() & 0xF;

        // initialize plaintexts (row 0)
        for (int j = 0; j < 32; ++j)
        {
            pt1[0][j] = rand() & 0xF;
            pt2[0][j] = pt1[0][j] ^ pt_delta[j];
        }

        enc_Step(pt1, ct1, key, num_rounds);
        enc_Step(pt2, ct2, key, num_rounds);
        for (size_t i = 0; i < 32; ++i)
        {
            ct3[num_rounds][i] = ct1[num_rounds][i] ^ ct_delta[i];
            ct4[num_rounds][i] = ct2[num_rounds][i] ^ ct_delta[i];
        }
        dec_Step(pt3, ct3, key, num_rounds);
        dec_Step(pt4, ct4, key, num_rounds);

        bool trial_successful = true;
        for (size_t i = 0; i < 32; ++i)
        {
            if ((ct3[0][i] ^ ct4[0][i]) != pt_delta[i])
            {
                trial_successful = false;
                break;
            }
        }
        successful_cnt += trial_successful ? 1 : 0;

        if (trial_successful)
        {
            for (int r = 0; r < num_rounds; r++) {
                if (check_all_equal_xor(pt1, pt2, ct3, ct4, r) == 0) {
                    check_cnt++;
                    cout << "find" << endl;
                    break;
                }
            }
        }

        delete[] pt1;
        delete[] pt2;
        delete[] pt3;
        delete[] pt4;
        delete[] ct1;
        delete[] ct2;
        delete[] ct3;
        delete[] ct4;

    }

    double prob1 = (double)successful_cnt / (double)nums;

    double prob2 = (double)check_cnt / (double)successful_cnt;


    //cout << p1 << ' ' << p2 << "  ";
    //write << p1 << ' ' << p2 << ": 2^" << log2(prob1) << ' ' << check_cnt << endl;
    cout << successful_cnt << ' ' << check_cnt << endl;
    printf("overall: 2^%.3f \n", log2(prob1));
    printf("overall: 2^%.3f \n", log2(prob2));
}