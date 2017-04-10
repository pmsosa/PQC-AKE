#ifndef KEM_H
#define KEM_H

#include "params.h"
//#include "FFT.h"
#include "Random.h"
#include "Algebra.h"
#include "KEM.h"


//void modCoeffs(ZZX& f, ZZ p);
void KEMKeyGen(ZZX& Kd, ZZX& Ke, ZZX& Kd_inv2);
void Encapsulate(ZZX& Ke, ZZX& c, ZZX& k);
void Decapsulate(ZZX& Kd, ZZX& c, ZZX& k, ZZX& Kd_inv2);
void run_KEM_example();
#endif