#ifndef DIGITALSIGNATURE_H
#define DIGITALSIGNATURE_H

#include "params.h"
#include "Sampling.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"


void SigKeyGen(ZZX Ks[2],ZZ_pX& Kv, MSK_Data* MSKD);
void Sign(ZZX s[2],vec_ZZ& msg, vec_ZZ& r, MSK_Data* MSKD);
bool Verify(ZZX Kv,ZZX s[2], vec_ZZ& msg, vec_ZZ& r);
void Hash(vec_ZZ&  hashed, vec_ZZ& r, vec_ZZ& msg);
// void modCoeffs(ZZX& f, ZZ p);

void Sign2(ZZX s[2],vec_ZZ& m2, vec_ZZ& msg, MSK_Data* MSKD);
bool Verify2(ZZX Kv,ZZX s[2], vec_ZZ& m2, vec_ZZ& m1);

void run_DS_example();
void run_DS_exampleMR();


void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey);
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey);
void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD);
void CompleteMSK(MSK_Data * MSKD, ZZX * MSK);
//void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK);
void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data * const MSKD);
unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id, const MSK_Data * const MSKD);
//void IBE_Encrypt(long C[2][N0], const long m[N0], const long id0[N0], const MPK_Data * const MPKD);
//void IBE_Decrypt(long message[N0], const long C[2][N0], const CC_t * const SKid_FFT);
//void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD);
//void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);
//void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD);
//void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);
//void DigSig2(const unsigned int nb_extr, MSK_Data * MSKD);
#endif
