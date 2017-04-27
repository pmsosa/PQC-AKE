#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>



//#include "Sampling.h"
#include "params.h"
//#include "FFT.h"
#include "Random.h"
#include "Algebra.h"
#include "KEM.h"
#include "FFT.h"

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();
const ZZ kem_norm = conv<ZZ>(KEM_NORM); //Move this to params


const bool dtime = false;  //Print Timing info?
const bool debug = false; //Print Debug info?



void KEMKeyGen(ZZX& Kd, ZZX& Ke, ZZX& Kd_inv2){
    

    ZZX f,g,g_inv,g_inv2,f_inv,f_inv2;
    bool fvalid = false;
    bool gvalid = false;
    int fretries = 0;
    int gretries = 0;

    while (!fvalid){

        try{
            //Build f & g        
            f = RandomPolyFixedSqNorm(kem_norm,N0-1);
            
            //Check that f and g are invertible in both Zq[x]/<x^n+1> and Z2[x]/<x^n+1>
            f_inv2 = conv<ZZX>(Inverse2(conv<ZZX>(f),2));     
            f_inv = conv<ZZX>(Inverse2(conv<ZZX>(f),q0));
 
            fvalid = true;
            
        }
        catch (int e){
            fretries += 1;
            //Keep going!
        }
    }


    while (!gvalid){

        try{
            //Build f & g        
            g = RandomPolyFixedSqNorm(kem_norm,N0-1);
            
            

            //Check that f and g are invertible in both Zq[x]/<x^n+1> and Z2[x]/<x^n+1>
            g_inv2 = conv<ZZX>(Inverse2(conv<ZZX>(g),2));            
            g_inv = conv<ZZX>(Inverse2(conv<ZZX>(g),q0));

             
            gvalid = true;
            
        }
        catch (int e){
            gretries += 1;
            //Keep going!
        }
    }


    //Debuggin' Information
    if (debug){
        cout << "\n\nKEM - Keygen\n";
        cout << "retries:" << fretries+gretries << "\n";
        cout << "f: "<< f << " ("<<deg(f)<<")"<<"\n";
        cout << "g: "<< g << " ("<<deg(g)<<")"<<"\n";

        cout << "f^-1: "<< f_inv << " ("<<deg(f_inv)<<")"<<"\n";
        cout << "f^-1 (2): "<< f_inv2 << " ("<<deg(f_inv2)<<")"<<"\n";

        cout << "g^-1: "<< g_inv << " ("<<deg(g_inv)<<")"<<"\n";
        cout << "g^-1 (2): "<< g_inv2 << " ("<<deg(g_inv2)<<")"<<"\n";
    }

    Kd = g;
    Kd_inv2 = g_inv2;
    Ke = (f*g_inv)%phi;
    //cout << "!"<<Ke<<"!";
    modCoeffs(Ke,q1);

    return;
}


void Encapsulate(ZZX& Ke, ZZX& c, ZZX& k){

    ZZ_p::init(conv<ZZ>(q0));

    ZZX r = RandomPolyFixedSqNorm(kem_norm,N0-1);
    ZZX e = RandomPolyFixedSqNorm(kem_norm,N0-1);

    FFTmultiply(c,Ke,r);
    c = conv<ZZX>((2*conv<ZZ_pX>(c))+conv<ZZ_pX>(e))%phi;

    // c = conv<ZZX>((2*conv<ZZ_pX>(Ke)*conv<ZZ_pX>(r))+conv<ZZ_pX>(e))%phi;

    modCoeffs(c,q1);

    ZZ_p::init(conv<ZZ>(2));
    k = conv<ZZX>(conv<ZZ_pX>(e));

    //Debuggin' Information
    if (debug){
        cout << "\n\nKEM - Enc\n";
        cout << "r: " << r << "\n";
        cout << "e: " << e << "\n"; 

        cout << "c: " << c << "\n";
        cout << "k: " << k << "\n"; 
    }
    return;
}

void Decapsulate(ZZX& Kd, ZZX& c, ZZX& k, ZZX& Kd_inv2){

    ZZ_p::init(conv<ZZ>(q0));

    FFTmultiply(k,Kd,c);
    k = conv<ZZX>(conv<ZZ_pX>(k))%phi;
    modCoeffs(k,q1);

    //cout << k;

    //ZZ_p::init(conv<ZZ>(2));

    //ZZX Kd_inv = conv<ZZX>(Inverse2(conv<ZZX>(Kd),2));
    modCoeffs(k,conv<ZZ>(2));
    //This can also be done by conv k to ZZ_pX
    
    //cout << "!!"<<k<<"!!\n";
    //FFT THIS!
    //cout << deg(k) << " -- " << deg(Kd_inv2) << "\n";
    FFTmultiply(k,k,Kd_inv2);
    k = (k)%phi;//conv<ZZX>(conv<ZZ_pX>(k)/conv<ZZ_pX>(Kd));
    modCoeffs(k,conv<ZZ>(2));

    return;
}


void run_KEM_example(){

    clock_t t1, t2;
    float t_keygen, t_enc, t_dec;

    ZZX Kd,Ke,Kd_inv2;
    
    t1 = clock();
    KEMKeyGen(Kd,Ke,Kd_inv2);
    t2 = clock();
    t_keygen = ((float)t2 - (float)t1)/1000000.0F;
 
    cout << "\nKd: " << Kd << " | Ke:" << Ke << "\n";

    t1 = clock();
    ZZX c,k;
    Encapsulate(Ke,c,k);
    t2 = clock();
    t_enc = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\nc: " << c << " | k:" << k << "\n";

    t1 = clock();
    ZZX k2;
    Decapsulate(Kd,c,k2,Kd_inv2);
    t2 = clock();
    t_dec = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\nk':" << k2 << "\n";

    bool valid = (k2 == k);
    cout << "VALID (k'==k): "<<valid << "\n";

    if (dtime){
        cout << "\nTiming\n";
        cout << "KEMKeyGen: " << t_keygen << "\n";
        cout << "KEMEnc   : " << t_enc << "\n";
        cout << "KEMDec   : " << t_dec << "\n";
    }


    return;
}


// clock_t t1, t2;
// t1 = clock();
//     t2 = clock();
//     diff = ((float)t2 - (float)t1)/1000000.0F;
//     cout << "\n\nIt took " << diff << " seconds to extract " << nb_extr << " keys." << endl;
//     cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per key." << endl << endl;
// }