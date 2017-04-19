#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>
#include <openssl/sha.h>

#include "Sampling.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"
#include "DigitalSignature.h"

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();


const bool dtime = false;  //Print Timing Info
const bool debug = false;  //Print Debug Info

///DIGITAL SIGNATURE - Pedro M. Sosa///

void SigKeyGen(ZZX Ks[2],ZZ_pX& Kv, MSK_Data* MSKD){

    
    ZZX MSK[4];
    ZZ_pX MPK;


    //GENERATE MATRIX B
    Keygen(MPK, MSK);
    CompleteMSK(MSKD, MSK);

    //ZZX f = MSKD -> PrK[0]; //Notice both paper invert f <==> g
    //ZZX g = MSKD -> PrK[1]; //So for simplicity think of them as nom/denom

    Ks[0] = MSKD -> PrK[0];   // denom
    Ks[1] = MSKD -> PrK[1]; // nom

    Kv = Quotient(Ks[0],Ks[1]); //h

    return;
}

/////////////////////////////////////////////////////////////////////
// Scheme with message recovery!
/////////////////////////////////////////////////////////////////////

void Sign(ZZX s[2],vec_ZZ& msg, vec_ZZ& r, MSK_Data* MSKD){
    //Hash Function
    r = conv<vec_ZZ>(RandomPoly(50));
    vec_ZZ hashed;
    Hash(hashed, r, msg);
    IBE_Extract(s, hashed, MSKD);

    return;
}

bool Verify(ZZX Kv,ZZX s[2], vec_ZZ& msg, vec_ZZ& r){
    //Hash id
    const ZZX phi = Cyclo();
    vec_ZZ hashed = vec_ZZ();
    Hash(hashed, r, msg);
    ZZX c;
    FFTmultiply(c,Kv,s[1]);
    ZZ_pX aux = conv<ZZ_pX>((conv<ZZX>(hashed) - (c))%phi);

    s[0] = conv<ZZX>(aux);
    modCoeffs(s[0],q1);

    double norm1 = 0;
    double norm2 = 0;
    double norm3 = 0;
    for(int i=0; i < deg(s[0]);i++){
        norm1 += conv<double>(s[0][i]*s[0][i]);
        norm2 += conv<double>(s[1][i]*s[1][i]);
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    norm3 = sqrt(norm1*norm1+norm2*norm2);
    if (debug){
        cout << norm3 << " |" << conv<ZZ>(DS_SIGMA*sqrt(q0/2*N0)*sqrt(2*N0));
    }
    return (norm3 < conv<ZZ>(DS_SIGMA*sqrt(q0/2*N0)*sqrt(2*N0)));//conv<ZZ>(1.36*q0/2*sqrt(2*N0)));
}

void Hash(vec_ZZ&  hashed, vec_ZZ& r, vec_ZZ& msg){
    SHA256_CTX ctx;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Init(&ctx);

    //vec_ZZ msg = RandomVector();
    hashed = msg;

    // if (debug){
    //     cout << "msg:" << msg << "\n";
    //     cout << "  r:" << r << "\n";
    //     cout << "m+r:" << hashed << "\n";
    // }

    for(int i = 0; i < msg.length(); i++) {
        char f = (conv<int>(msg[i])%255);
        SHA256_Update(&ctx, &f, 1);
    }
    for(int i = 0; i < r.length(); i++) {
        char f = (conv<int>(r[i])%255);
        SHA256_Update(&ctx, &f, 1);
    }

    SHA256_Final(digest, &ctx);


    for(int i=0; i < SHA256_DIGEST_LENGTH and i < N0; i++){
        //cout << int(digest[i]) << " ";
        hashed[i] = conv<ZZ>(digest[i])%q0;
    }
    if (debug){
        cout << "Hashed: " << hashed << "\n";
    }
}

/////////////////////////////////////////////////////////////////////
// Scheme without message recovery!
/////////////////////////////////////////////////////////////////////


void HashF(vec_ZZ&  hashed, vec_ZZ& msg){
    SHA256_CTX ctx;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Init(&ctx);
    //SHA256_Update(&ctx, "F", 1);

    //vec_ZZ msg = RandomVector();
    hashed.SetLength(32);
    //cout << SHA256_DIGEST_LENGTH <<"\n";
    //cout << "hashed:  ------ "<< hashed <<"\n\n";
    //cout << "MSG" << msg << "\n";

    // if (debug){
    //     cout << "msg:" << msg << "\n";
    //     cout << "  r:" << r << "\n";
    //     cout << "m+r:" << hashed << "\n";
    // }

    //cout << "Elem: ";
    for(int i = 0; i < msg.length(); i++) {
        char f = (conv<int>(msg[i])%255);
        SHA256_Update(&ctx, &f, 1);
        //cout << int(f) << ",";
    }

    //cout << "\nHASH:";
    SHA256_Final(digest, &ctx);
    for(int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        //cout << int(digest[i]) <<",";
    }


    for(int i=0; i < 32 and i < N0; i++){
        //cout << int(digest[i]) << " ";
        hashed[i] = conv<ZZ>(digest[i])%q0;
    }
    //cout << "\n\nh:" <<hashed <<"\n";
    if (debug){
        cout << "Hashed F: " << hashed << "\n";
    }
}


void HashH(vec_ZZ&  hashed, vec_ZZ& msg){
    SHA256_CTX ctx;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Init(&ctx);
    SHA256_Update(&ctx, "H", 1);

    //vec_ZZ msg = RandomVector();
    hashed.SetLength(32);

    // if (debug){
    //     cout << "msg:" << msg << "\n";
    //     cout << "  r:" << r << "\n";
    //     cout << "m+r:" << hashed << "\n";
    // }

    for(int i = 0; i < msg.length(); i++) {
        char f = (conv<int>(msg[i])%255);
        SHA256_Update(&ctx, &f, 1);
    }

    SHA256_Final(digest, &ctx);


    for(int i=0; i < 32 and i < N0; i++){
        //cout << int(digest[i]) << " ";
        hashed[i] = conv<ZZ>(digest[i])%q0;
    }
    if (debug){
        cout << "Hashed H: " << hashed << "\n";
    }
}

void Sign2(ZZX s[2],vec_ZZ& m2, vec_ZZ& msg, MSK_Data* MSKD){
    vec_ZZ m1;
    vec_ZZ h1,f1,t;

    // 512: |m1|=496, 1024: |m1|= 1008

    int m1_size = 480;
    m1.SetLength(m1_size);
    m2.SetLength(msg.length()-m1_size);

    for (int i=0; i < m1_size; i++){
        m1[i] = msg[i];
        if (i+m1_size < msg.length()){
            m2[i] = msg[i+m1_size];
        }
    }
    if (debug){
        cout << "M:"  << msg << "\n";
        cout << "M1:" << m1  << "\n";
        cout << "M2:" << m2  << "\n";
    }

    //cout << "M1:" << m1  << "\n";

    // t = (m1 + F(H(m))) mod q || H(m))
    HashH(h1,msg);
    HashF(f1,h1);
    //cout << "f1:"<<f1 <<"\n";
    //cout << "h1:"<<h1 <<"\n";
    //cout << "f1:"<<f1 << "\n";
    f1.SetLength(480);
    for (int i=0; i < f1.length(); i++){
        f1[i] = (f1[i]+m1[i])%q0;
    }



    cout << "To Be Recovered:"<<m1 <<"\n";

    t = f1;
    t.append(h1);

    //cout << "t:"<<t <<"\n";
    //cout << "\nt:"<< t.length()<<"\n";


   

    //cout << "\nt:"<< t.length()<<"\n";
    
    // s1,s2 DS such that hs1 +s2 = t mod q
    IBE_Extract(s, t, MSKD);





    // This returns s[1] and s[2] and m2 (last part of message)
}

bool Verify2(ZZX Kv,ZZX s[2], vec_ZZ& m2, vec_ZZ& m1){
    //Hash id
    vec_ZZ m,t,t1,h1,f1,t2;

    const ZZX phi = Cyclo();
    vec_ZZ hashed = vec_ZZ();
    ZZX c;
    FFTmultiply(c,Kv,s[1]);
    c = c + s[0];
    t = conv<vec_ZZ>(c);

    t1.SetLength(480);
    h1.SetLength(32);



    for (int i=0; i < 480; i++){
        t1[i] = t[i]%q0;
        if (i < 32){
            h1[i] = t[i+480]%q0;
        }

    }



    //cout << "t:" <<t  <<"\n";
    //cout << "t1:"<<t1 <<"\n";
    //cout << "h1:"<<h1<<"\n";
    //cout << "m2o:"<<m2_o<<"\n";

    HashF(f1,h1);

    //cout << "f1:"<<f1 <<"\n";

    m1.SetLength(t1.length());
    f1.SetLength(t1.length());
    for (int i=0; i < t1.length(); i++){
        m1[i] = (t1[i] - f1[i])%q0;
    }

    m = m1;
    m.append(m2);
    //cout << "m: " << m1 << "\n";
    HashH(t2,m1);

    //cout << "t2"<<t2 <<"\n";
    //cout << "h1"<<h1 <<"\n";
    for (int i=0; i < h1.length(); i++){
        if (h1[i] != t2[i]){
            return false;
        }
    }
    // m1 = conv<vec_ZZ>(t);
    // modCoeffs(s[0],q1);

    double norm1 = 0;
    double norm2 = 0;
    double norm3 = 0;
    for(int i=0; i < deg(s[0]);i++){
        norm1 += conv<double>(s[0][i]*s[0][i]);
        norm2 += conv<double>(s[1][i]*s[1][i]);
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    norm3 = sqrt(norm1*norm1+norm2*norm2);
    if (debug){
        cout << norm3 << " |" << conv<ZZ>(DS_SIGMA*sqrt(q0/2*N0)*sqrt(2*N0));
    }
    return (norm3 < conv<ZZ>(DS_SIGMA*sqrt(q0/2*N0)*sqrt(2*N0)));//conv<ZZ>(1.36*q0/2*sqrt(2*N0)));
}





void run_DS_example(){


    clock_t t1, t2;
    float t_keygen, t_sig, t_ver;

    /// RUNNING THE KEYGEN ALGORITHM
        
        ZZX Ks[2];
        ZZ_pX Kv;
        MSK_Data * MSKD = new MSK_Data;

        t1 = clock();
        SigKeyGen(Ks,Kv,MSKD);
        t2 = clock();
        t_keygen = ((float)t2 - (float)t1)/1000000.0F;
        //ZZX Kv2 = conv<ZZX>(Kv);

    /// RUNNING THE SIGN ALGORITHM
        vec_ZZ r,m2;
        vec_ZZ msg = RandomVector();
        ZZX s[2];
        t1 = clock();
        Sign2(s,m2,msg,MSKD);

        t2 = clock();
        t_sig = ((float)t2 - (float)t1)/1000000.0F;

        //PRINTIN DEBUG INFO
        if (debug){
            cout << "\n---Sign---\n";
            cout << "\n (rand) message: ";
            for(int i = 0; i < N0;i++){
                cout << msg[i] << ",";
            }
            cout <<"\n";

            cout << "             s1: ";
            for(int i = 0; i < N0;i++){
                cout << s[0][i] << ",";
            }
            cout <<"\n";

            cout << "     (debug) s2: ";
            for(int i = 0; i < N0;i++){
                cout << s[1][i] << ",";
            }
            cout <<"\n";
        }

    /// RUNNING THE VERIFY ALGORITHM
        

        //Notice
        //1. Turn Kv into ZZX(Kv)
        //2. We will pretend that s[0] is empty since we didn't get that info
        vec_ZZ recovery;
        t1 = clock();
        bool valid = Verify2(conv<ZZX>(Kv),s, m2,recovery);
        t2 = clock();
        t_ver = ((float)t2 - (float)t1)/1000000.0F;

        cout <<"Recovery" << recovery;

        //Not Certain of this part.
        // double norm1 = 0;
        // double norm2 = 0;
        // double norm3 = 0;
        // for(int i=0; i < deg(s[0]);i++){
        //     norm1 += conv<double>(s[0][i]*s[0][i]);
        //     norm2 += conv<double>(s[1][i]*s[1][i]);
        // }
        // norm1 = sqrt(norm1);
        // norm2 = sqrt(norm2);
        // norm3 = sqrt(norm1*norm1+norm2*norm2);
        // cout << norm3 << " | " << conv<ZZ>(1.36*q0/2*sqrt(2*N0));

        //PRINTIN DEBUG INFO
        if (debug){
            cout << "\n--Verify---\n";
            cout << "     derived s1: ";
            for (int i=0; i < N0; i++){
                cout << s[0][i] << ",";
            }
            cout <<"\n";

            cout << "          Valid: " << valid <<"\n" ;
        }

    if (dtime){
        cout << "\nTiming\n";
        cout << "DSKeyGen : " << t_keygen << "\n";
        cout << "DSSign   : " << t_sig << "\n";
        cout << "DSVerif  : " << t_ver << "\n";
    }

    return;
}


void run_DS_example_orig(){


    clock_t t1, t2;
    float t_keygen, t_sig, t_ver;

    /// RUNNING THE KEYGEN ALGORITHM
        
        ZZX Ks[2];
        ZZ_pX Kv;
        MSK_Data * MSKD = new MSK_Data;

        t1 = clock();
        SigKeyGen(Ks,Kv,MSKD);
        t2 = clock();
        t_keygen = ((float)t2 - (float)t1)/1000000.0F;
        //ZZX Kv2 = conv<ZZX>(Kv);

        //PRINTIN DEBUG INFO
        if (debug){
            cout << "\n---Keygen---\n";
            cout << " Private:\n";
            cout << "    denom: ";
            for( int i = 0; i< N0;i++){
                cout << Ks[0][i] << ",";
            }
            cout <<"\n";

            cout <<"      nom: ";
            for (int i = 0; i< N0;i++){
                cout << Ks[1][i] << ",";
            }
            cout <<"\n";

            cout << "\n Public:\n";

            cout << "        h: ";
            for(int i = 0; i < N0;i++){
                cout << Kv[i] << ","; //"( "<< Kv2[i] <<" ) " << ",";
            }
            cout <<"\n";
        }
    /// RUNNING THE SIGN ALGORITHM
        vec_ZZ r;
        vec_ZZ msg = RandomVector();
        ZZX s[2];
        t1 = clock();
        Sign(s,msg,r,MSKD);
        t2 = clock();
        t_sig = ((float)t2 - (float)t1)/1000000.0F;

        //PRINTIN DEBUG INFO
        if (debug){
            cout << "\n---Sign---\n";
            cout << "\n (rand) message: ";
            for(int i = 0; i < N0;i++){
                cout << msg[i] << ",";
            }
            cout <<"\n";

            cout << "             s1: ";
            for(int i = 0; i < N0;i++){
                cout << s[0][i] << ",";
            }
            cout <<"\n";

            cout << "     (debug) s2: ";
            for(int i = 0; i < N0;i++){
                cout << s[1][i] << ",";
            }
            cout <<"\n";
        }

    /// RUNNING THE VERIFY ALGORITHM
        

        //Notice
        //1. Turn Kv into ZZX(Kv)
        //2. We will pretend that s[0] is empty since we didn't get that info
        t1 = clock();
        bool valid = Verify(conv<ZZX>(Kv),s, msg,r);
        t2 = clock();
        t_ver = ((float)t2 - (float)t1)/1000000.0F;


        //Not Certain of this part.
        // double norm1 = 0;
        // double norm2 = 0;
        // double norm3 = 0;
        // for(int i=0; i < deg(s[0]);i++){
        //     norm1 += conv<double>(s[0][i]*s[0][i]);
        //     norm2 += conv<double>(s[1][i]*s[1][i]);
        // }
        // norm1 = sqrt(norm1);
        // norm2 = sqrt(norm2);
        // norm3 = sqrt(norm1*norm1+norm2*norm2);
        // cout << norm3 << " | " << conv<ZZ>(1.36*q0/2*sqrt(2*N0));

        //PRINTIN DEBUG INFO
        if (debug){
            cout << "\n--Verify---\n";
            cout << "     derived s1: ";
            for (int i=0; i < N0; i++){
                cout << s[0][i] << ",";
            }
            cout <<"\n";

            cout << "          Valid: " << valid <<"\n" ;
        }

    if (dtime){
        cout << "\nTiming\n";
        cout << "DSKeyGen : " << t_keygen << "\n";
        cout << "DSSign   : " << t_sig << "\n";
        cout << "DSVerif  : " << t_ver << "\n";
    }

    return;
}



///Repurposed IBE Code - Peikart///


//==============================================================================
//Generates from parameters N and q :
// - a public key : polynomial h
// - a private key : polynomials f,g,F,G
//==============================================================================
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey)
{
    ZZ SqNorm;
    ZZX f,g,F,G;

    SqNorm = conv<ZZ>(1.36*q0/2);

    GenerateBasis(f, g, F, G, SqNorm);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for(unsigned int i=0; i<4; i++)
    {
            PrivateKey[i].SetLength(N0);
    }

    PublicKey = Quotient(f, g);
}

//==============================================================================
//Computes the private basis B from private key PrivateKey and parameter N
//==============================================================================
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey)
{
    ZZX f,g,F,G;
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;

    B = BasisFromPolynomials(g, f, G, F);
}


void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD)
{

    int i;
    unsigned j;
    RR_t ci[2*N0], zi, cip, sip, aux;

    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    }

    for(j=0; j<2*N0; j++)
    {

    }    

    for(i=2*N0-1; i>=0; i--)
    {
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, MSKD->Bstar[i])/(aux*aux);
        sip = s/aux;
        zi = Sample4(cip, sip*PiPrime);

        for(j=0; j<2*N0; j++)
        {
            ci[j] -= zi*(MSKD->B)[i][j];
        }
    }

    for(j=0; j<2*N0; j++)
    {
        v[j] = c[j] - ci[j];
    }

}



//==============================================================================
//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================
//==============================================================================


void CompleteMSK(MSK_Data * MSKD, ZZX * MSK)
{
    unsigned int i, j;
    mat_ZZ B0;

    for(i=0; i<4; i++)
    {
        MSKD->PrK[i] = MSK[i];
        ZZXToFFT(MSKD->PrK_fft[i], MSK[i]);
    }

    CompletePrivateKey(B0, MSK);

    for(i=0; i<2*N0; i++)
    {
        for(j=0; j<2*N0; j++)
        {
            MSKD->B[i][j] = ( (RR_t) conv<double>(B0[i][j]) );
        }
    }

    for(i=0; i<1; i++)
    {
        FastMGS(MSKD->Bstar, MSKD->B);
    }

    for(i=0; i<2*N0; i++)
    {
        MSKD->GS_Norms[i] = sqrt( DotProduct(MSKD->Bstar[i], MSKD->Bstar[i]) );
    }

    MSKD->sigma = 2*MSKD->GS_Norms[0];

}

void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data * const MSKD)
{
    unsigned int i;
    RR_t c[2*N0], sk[2*N0], sigma;
    ZZX f,g,aux;

    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    sigma = MSKD->sigma;
    SK_id[0].SetLength(N0);
    SK_id[1].SetLength(N0);

    for(i=0;i<N0;i++)
    {
        c[i] = ((RR_t) conv<double>(id[i])) ;
        c[i+N0] = 0;
    }

    GPV(sk, c, sigma, MSKD);

    for(i=0; i<N0; i++)
    {
        sk[i] = c[i] - sk[i];
        sk[i+N0] = - sk[i+N0];
    }

    for(i=0; i<N0; i++)
    {
        SK_id[0][i] = sk[i];
        SK_id[1][i] = sk[i+N0];
    }
    
}