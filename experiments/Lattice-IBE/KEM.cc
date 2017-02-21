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

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();
const ZZ kem_norm = conv<ZZ>(30); //Move this to params


void modCoeffs(ZZX& f, ZZ p){
    ZZ pp = p/2;
    for (int i=0;i <= deg(f);i++){
        ZZ temp = f[i]%p;
        if (temp>pp){ temp = temp - p;}

        SetCoeff(f,i,temp);
    }
}

void KEMKeyGen(ZZX& Kd, ZZX& Ke){
    

    bool valid = false;
    int retries = 0;
    while (!valid){

        try{
            //Build f & g        
            ZZX f = RandomPolyFixedSqNorm(kem_norm,N0-1);
            ZZX g = RandomPolyFixedSqNorm(kem_norm,N0-1);
            
            

            //Check that f is invertible in both Zq[x]/<x^n+1> and Z2[x]/<x^n+1>
            ZZX f_inv = conv<ZZX>(Inverse2(conv<ZZX>(f),q0));
            ZZX f_inv2 = conv<ZZX>(Inverse2(conv<ZZX>(f),2));

            //Check that g is invertible in both Zq[x]/<x^n+1> and Z2[x]/<x^n+1>
            ZZX g_inv = conv<ZZX>(Inverse2(conv<ZZX>(g),q0));
            ZZX g_inv2 = conv<ZZX>(Inverse2(conv<ZZX>(g),2));


            //Debuggin' Information
            cout << "\n\nKEM - Keygen\n";
            cout << "retries:" << retries << "\n";
            cout << "f: "<< f << " ("<<deg(f)<<")"<<"\n";
            cout << "g: "<< g << " ("<<deg(g)<<")"<<"\n";

            cout << "f^-1: "<< f_inv << " ("<<deg(f_inv)<<")"<<"\n";
            cout << "f^-1 (2): "<< f_inv2 << " ("<<deg(f_inv2)<<")"<<"\n";

            cout << "g^-1: "<< g_inv << " ("<<deg(g_inv)<<")"<<"\n";
            cout << "g^-1 (2): "<< g_inv2 << " ("<<deg(g_inv2)<<")"<<"\n";
            valid = true;


            Kd = g;
            Ke = (f*g_inv)%phi;
            //cout << "!"<<Ke<<"!";
            modCoeffs(Ke,q1);
            
        }
        catch (int e){
            retries += 1;
            //Keep going!
        }
    }

    return;
}


void Encapsulate(ZZX& Ke, ZZX& c, ZZX& k){

    ZZ_p::init(conv<ZZ>(q0));

    ZZX r = RandomPolyFixedSqNorm(kem_norm,N0-1);
    ZZX e = RandomPolyFixedSqNorm(kem_norm,N0-1);

    c = conv<ZZX>((2*conv<ZZ_pX>(Ke)*conv<ZZ_pX>(r))+conv<ZZ_pX>(e))%phi;

    modCoeffs(c,q1);

    ZZ_p::init(conv<ZZ>(2));
    k = conv<ZZX>(conv<ZZ_pX>(e));

    //Debuggin' Information
    cout << "\n\nKEM - Enc\n";
    cout << "r: " << r << "\n";
    cout << "e: " << e << "\n"; 

    cout << "c: " << c << "\n";
    cout << "k: " << k << "\n"; 

    return;
}

void Decapsulate(ZZX& Kd, ZZX& c, ZZX& k){

    ZZ_p::init(conv<ZZ>(q0));

    k = conv<ZZX>(conv<ZZ_pX>(Kd*c))%phi;
    modCoeffs(k,q1);

    //cout << k;

    //ZZ_p::init(conv<ZZ>(2));

    ZZX Kd_inv = conv<ZZX>(Inverse2(conv<ZZX>(Kd),2));
    modCoeffs(k,conv<ZZ>(2));
    //This can also be done by conv k to ZZ_pX
    
    //cout << "!!"<<k<<"!!\n";
    
    k = (k*Kd_inv)%phi;//conv<ZZX>(conv<ZZ_pX>(k)/conv<ZZ_pX>(Kd));
    modCoeffs(k,conv<ZZ>(2));

    return;
}


void run_KEM_example(){

    ZZX Kd,Ke;
    KEMKeyGen(Kd,Ke);
    cout << "\nKd: " << Kd << " | Ke:" << Ke << "\n";

    ZZX c,k;
    Encapsulate(Ke,c,k);
    cout << "\nc: " << c << " | k:" << k << "\n";


    ZZX k2;
    Decapsulate(Kd,c,k2);
    cout << "\nk':" << k2 << "\n";

    bool valid = (k2 == k);
    cout << "VALID (k'==k): "<<valid << "\n";


    return;
}