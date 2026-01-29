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

typedef uint16_t u16;
typedef uint8_t u8;


#define ROUND0 20
#define ROUND1 8

static const u8 S[16] = { 0x0E, 0x09, 0x0F, 0x00, 0x0D, 0x04, 0x0A, 0x0B, 0x01, 0x02, 0x08, 0x03, 0x07, 0x06, 0x0C, 0x05 };

void KeySchedule(int nrounds, u8 key[10], u8 output[][4])
{
    u8 i, KeyR[3];

    output[0][3] = key[9];
    output[0][2] = key[8];
    output[0][1] = key[7];
    output[0][0] = key[6];

    for (i = 1; i < nrounds; i++)
    {
        // K <<< 24            
        KeyR[0] = key[7];
        KeyR[1] = key[8];
        KeyR[2] = key[9];
        // 
        key[9] = key[6];
        key[8] = key[5];
        key[7] = key[4];
        key[6] = key[3];
        key[5] = key[2];
        key[4] = key[1];
        key[3] = key[0];
        // 
        key[2] = KeyR[2];
        key[1] = KeyR[1];
        key[0] = KeyR[0];


        // reste du keyschedule                 
        key[6] ^= (S[(key[9] >> 4) & 0x0F] << 4);
        key[3] ^= (S[key[9] & 0x0F] << 4);
        key[8] ^= ((key[8] >> 4) & 0x0F);
        key[6] ^= ((key[1] >> 4) & 0x0F);

        key[6] ^= ((i & 0x1F) << 2);

        output[i][3] = key[9];
        output[i][2] = key[8];
        output[i][1] = key[7];
        output[i][0] = key[6];
    }
}

void Swap(u8 block[8])
{
    u8 tmp[4];

    tmp[0] = block[0];
    tmp[1] = block[1];
    tmp[2] = block[2];
    tmp[3] = block[3];

    block[0] = block[4];
    block[1] = block[5];
    block[2] = block[6];
    block[3] = block[7];

    block[4] = tmp[0];
    block[5] = tmp[1];
    block[6] = tmp[2];
    block[7] = tmp[3];
}

void OneRound_Step(u8 x[8], u8 k[4], u8 y[8])
{
    u8 t[4], tmp[4];

    // AJOUT CLE
    tmp[0] = x[4] ^ k[0];
    tmp[1] = x[5] ^ k[1];
    tmp[2] = x[6] ^ k[2];
    tmp[3] = x[7] ^ k[3];

    // PASSAGE DANS LES BOITES S
    tmp[0] = ((S[((tmp[0]) >> 4) & 0x0F]) << 4) ^ S[(tmp[0] & 0x0F)];
    tmp[1] = ((S[((tmp[1]) >> 4) & 0x0F]) << 4) ^ S[(tmp[1] & 0x0F)];
    tmp[2] = ((S[((tmp[2]) >> 4) & 0x0F]) << 4) ^ S[(tmp[2] & 0x0F)];
    tmp[3] = ((S[((tmp[3]) >> 4) & 0x0F]) << 4) ^ S[(tmp[3] & 0x0F)];

    // PASSAGE DE LA PERMUTATION P
    t[0] = ((tmp[0] >> 4) & 0x0F) ^ (tmp[1] & 0xF0);
    t[1] = (tmp[0] & 0x0F) ^ ((tmp[1] & 0x0F) << 4);
    t[2] = ((tmp[2] >> 4) & 0x0F) ^ (tmp[3] & 0xF0);
    t[3] = (tmp[2] & 0x0F) ^ ((tmp[3] & 0x0F) << 4);
    // FIN DE LA FONCTION F

    // PARTIE GAUCHE AVEC DECALAGE DE 8 SUR LA GAUCHE  
    tmp[0] = x[3] ^ t[0];
    tmp[1] = x[0] ^ t[1];
    tmp[2] = x[1] ^ t[2];
    tmp[3] = x[2] ^ t[3];

    // PARTIE DROITE
    y[0] = tmp[0];
    y[1] = tmp[1];
    y[2] = tmp[2];
    y[3] = tmp[3];
    y[4] = x[4];
    y[5] = x[5];
    y[6] = x[6];
    y[7] = x[7];
}

void OneRound(u8 x[8], u8 k[4])
{
    u8 t[4], tmp[4];

    // AJOUT CLE
    tmp[0] = x[4] ^ k[0];
    tmp[1] = x[5] ^ k[1];
    tmp[2] = x[6] ^ k[2];
    tmp[3] = x[7] ^ k[3];

    // PASSAGE DANS LES BOITES S
    tmp[0] = ((S[((tmp[0]) >> 4) & 0x0F]) << 4) ^ S[(tmp[0] & 0x0F)];
    tmp[1] = ((S[((tmp[1]) >> 4) & 0x0F]) << 4) ^ S[(tmp[1] & 0x0F)];
    tmp[2] = ((S[((tmp[2]) >> 4) & 0x0F]) << 4) ^ S[(tmp[2] & 0x0F)];
    tmp[3] = ((S[((tmp[3]) >> 4) & 0x0F]) << 4) ^ S[(tmp[3] & 0x0F)];

    // PASSAGE DE LA PERMUTATION P
    t[0] = ((tmp[0] >> 4) & 0x0F) ^ (tmp[1] & 0xF0);
    t[1] = (tmp[0] & 0x0F) ^ ((tmp[1] & 0x0F) << 4);
    t[2] = ((tmp[2] >> 4) & 0x0F) ^ (tmp[3] & 0xF0);
    t[3] = (tmp[2] & 0x0F) ^ ((tmp[3] & 0x0F) << 4);
    // FIN DE LA FONCTION F

    // PARTIE GAUCHE AVEC DECALAGE DE 8 SUR LA GAUCHE  
    tmp[0] = x[3] ^ t[0];
    tmp[1] = x[0] ^ t[1];
    tmp[2] = x[1] ^ t[2];
    tmp[3] = x[2] ^ t[3];

    // PARTIE DROITE
    x[0] = tmp[0];
    x[1] = tmp[1];
    x[2] = tmp[2];
    x[3] = tmp[3];

}

void Encrypt(int nrounds, u8 x[8], u8 y[8], u8 subkey[][4])
{
    u8 tmp[8];
    for (int i = 0; i < 8; i++) tmp[i] = x[i];

    for (int i = 0; i < nrounds; i++)
    {
        OneRound(tmp, subkey[i]);
        Swap(tmp);
    }

    for (int i = 0; i < 8; i++) y[i] = tmp[i];
}

void OneRound_Inv_Step(u8 y[8], u8 k[4], u8 x[8])
{
    u8 t[4], tmp[4];

    // FAIRE PASSER Y_0, Y_1, Y_2, Y_3 dans F
    // AJOUT CLE
    tmp[0] = y[4] ^ k[0];
    tmp[1] = y[5] ^ k[1];
    tmp[2] = y[6] ^ k[2];
    tmp[3] = y[7] ^ k[3];


    // PASSAGE DANS LES BOITES S
    tmp[0] = ((S[((tmp[0]) >> 4) & 0x0F]) << 4) ^ S[(tmp[0] & 0x0F)];
    tmp[1] = ((S[((tmp[1]) >> 4) & 0x0F]) << 4) ^ S[(tmp[1] & 0x0F)];
    tmp[2] = ((S[((tmp[2]) >> 4) & 0x0F]) << 4) ^ S[(tmp[2] & 0x0F)];
    tmp[3] = ((S[((tmp[3]) >> 4) & 0x0F]) << 4) ^ S[(tmp[3] & 0x0F)];

    // PASSAGE DE LA PERMUTATION P
    t[0] = ((tmp[0] >> 4) & 0x0F) ^ (tmp[1] & 0xF0);
    t[1] = (tmp[0] & 0x0F) ^ ((tmp[1] & 0x0F) << 4);
    t[2] = ((tmp[2] >> 4) & 0x0F) ^ (tmp[3] & 0xF0);
    t[3] = (tmp[2] & 0x0F) ^ ((tmp[3] & 0x0F) << 4);
    // FIN DE LA FONCTION F

    // PARTIE DROITE AVEC DECALAGE DE 8 SUR LA DROITE
    tmp[0] = y[0] ^ t[0];
    tmp[1] = y[1] ^ t[1];
    tmp[2] = y[2] ^ t[2];
    tmp[3] = y[3] ^ t[3];

    // PARTIE GAUCHE
    x[0] = tmp[1];
    x[1] = tmp[2];
    x[2] = tmp[3];
    x[3] = tmp[0];
    x[4] = y[4];
    x[5] = y[5];
    x[6] = y[6];
    x[7] = y[7];
}

void OneRound_Inv(u8 y[8], u8 k[4])
{
    u8 t[4], tmp[4];

    // FAIRE PASSER Y_0, Y_1, Y_2, Y_3 dans F
    // AJOUT CLE
    tmp[0] = y[4] ^ k[0];
    tmp[1] = y[5] ^ k[1];
    tmp[2] = y[6] ^ k[2];
    tmp[3] = y[7] ^ k[3];


    // PASSAGE DANS LES BOITES S
    tmp[0] = ((S[((tmp[0]) >> 4) & 0x0F]) << 4) ^ S[(tmp[0] & 0x0F)];
    tmp[1] = ((S[((tmp[1]) >> 4) & 0x0F]) << 4) ^ S[(tmp[1] & 0x0F)];
    tmp[2] = ((S[((tmp[2]) >> 4) & 0x0F]) << 4) ^ S[(tmp[2] & 0x0F)];
    tmp[3] = ((S[((tmp[3]) >> 4) & 0x0F]) << 4) ^ S[(tmp[3] & 0x0F)];

    // PASSAGE DE LA PERMUTATION P
    t[0] = ((tmp[0] >> 4) & 0x0F) ^ (tmp[1] & 0xF0);
    t[1] = (tmp[0] & 0x0F) ^ ((tmp[1] & 0x0F) << 4);
    t[2] = ((tmp[2] >> 4) & 0x0F) ^ (tmp[3] & 0xF0);
    t[3] = (tmp[2] & 0x0F) ^ ((tmp[3] & 0x0F) << 4);
    // FIN DE LA FONCTION F

    // PARTIE DROITE AVEC DECALAGE DE 8 SUR LA DROITE
    tmp[0] = y[0] ^ t[0];
    tmp[1] = y[1] ^ t[1];
    tmp[2] = y[2] ^ t[2];
    tmp[3] = y[3] ^ t[3];

    // PARTIE GAUCHE
    y[0] = tmp[1];
    y[1] = tmp[2];
    y[2] = tmp[3];
    y[3] = tmp[0];
}

void Decrypt(int nrounds, u8 y[8], u8 x[8], u8 subkey[][4])
{
    u8 tmp[8];
    for (int i = 0; i < 8; i++) tmp[i] = y[i];

    for (int i = nrounds - 1; i >= 0; i--)
    {
        Swap(tmp);
        OneRound_Inv(tmp, subkey[i]);
    }

    for (int i = 0; i < 8; i++) x[i] = tmp[i];
}

void Print4Bytes(u8 x[4]) {
    for (int i = 3; i >= 0; i--) {
        cout << hex << setw(2) << setfill('0') << static_cast<unsigned>(x[i]) << ' ';
    }
    cout << endl;
}

void Print8Bytes(u8 x[8]) {
    for (int i = 7; i >= 0; i--) {
        cout << hex << setw(2) << setfill('0') << static_cast<unsigned>(x[i]) << ' ';
    }
    cout << endl;
}

void Print10Bytes(u8 x[10]) {
    for (int i = 9; i >= 0; i--) {
        cout << hex << setw(2) << setfill('0') << static_cast<unsigned>(x[i]) << ' ';
    }
    cout << endl;
}

// Key to array
void hex_to_u8_array(char* s, u8* delta) {
    for (int i = 0; i < 8; i++) {
        char high = s[2 * i];
        char low = s[2 * i + 1];
        delta[7 - i] = (uint8_t)((strtol(string(1, high).c_str(), nullptr, 16) << 4) |
            strtol(string(1, low).c_str(), nullptr, 16));
    }
}

int main() {
    char pt_delta_str[] = "0000000000000000";
    char ct_delta_str[] = "0000000000000000"; 
    int num_rounds = 9; 

    pt_delta_str[7] = '4';
    //pt_delta_str[5] = '4';
    //ct_delta_str[4] = '4';
    ct_delta_str[8] = '4';

    u8 pt_delta[8], ct_delta[8];
    hex_to_u8_array(pt_delta_str, pt_delta);
    hex_to_u8_array(ct_delta_str, ct_delta);

    uint64_t nums_key = 1ULL << 6; // Num of Key
    uint64_t nums_msg = 1ULL << 36; // Num of messages per key 
   
    printf("LBlocks Boomerang experiment: 2^%.1f messages for each of 2^%.1f keys, rounds = %d\n",
        log2((double)nums_msg), log2((double)nums_key), num_rounds);
    printf("%s\n", pt_delta_str);
    printf("%s\n", ct_delta_str);

    uint64_t successful_cnt = 0;

#pragma omp parallel reduction(+:successful_cnt)
{
    std::mt19937 rng((unsigned)time(NULL) ^ omp_get_thread_num());
    std::uniform_int_distribution<int> dist(0, 255);

    u8 pt1[8], pt2[8], pt3[8], pt4[8];
    u8 ct1[8], ct2[8], ct3[8], ct4[8];

    u8 key[10], rkey[32][4];
    for (uint64_t ctr_key = 0; ctr_key < nums_key; ctr_key++) {
        for (int i = 0; i < 10; i++) key[i] = dist(rng);
        KeySchedule(num_rounds, key, rkey);

#pragma omp for schedule(dynamic,10000)

        for (uint64_t ctr_msg = 0; ctr_msg < nums_msg; ctr_msg++) {

            for (int i = 0; i < 8; i++) {
                pt1[i] = dist(rng);
                pt2[i] = pt1[i] ^ pt_delta[i];
            }

            Encrypt(num_rounds, pt1, ct1, rkey);
            Encrypt(num_rounds, pt2, ct2, rkey);

            for (int i = 0; i < 8; i++)
            {
                ct3[i] = ct1[i] ^ ct_delta[i];
                ct4[i] = ct2[i] ^ ct_delta[i];
            }

            Decrypt(num_rounds, ct3, pt3, rkey);
            Decrypt(num_rounds, ct4, pt4, rkey);

            bool ok = true;
            for (int i = 0; i < 8; i++) {
                if ((pt3[i] ^ pt4[i]) != pt_delta[i]) {
                    ok = false;
                    break;
                }
            }

            successful_cnt += ok;
        }
    }
}
    if (successful_cnt > 0) {
        double prob = (double)successful_cnt / (double)(nums_key * nums_msg);
        printf("Overall Probability = 2^%.3f\n", log2(prob));
    }
    else printf("No hits found.\n");

    return 0;
}