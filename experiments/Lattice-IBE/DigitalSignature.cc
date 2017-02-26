#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "Sampling.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"
#include "DigitalSignature.h"

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();


const bool dtime = true;  //Print Timing info?
const bool debug = true; //Print Debug info?

///DIGITAL SIGNATURE - Pedro M. Sosa///

// void InitDS(){

//     srand(rdtsc());                 //Seed Random

//     ZZ_p::init(q1);                 //q1 = ZZ(q0)

//     const ZZX phi = Cyclo();        //Phi - Used all over the place
//     ZZ_pX phiq = conv<ZZ_pX>(phi);  
//     ZZ_pXModulus PHI(phiq);
// }

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

void Sign(ZZX s[2],vec_ZZ msg, MSK_Data* MSKD){
    //TODO: HASH FUNCTION
    IBE_Extract(s, msg, MSKD);
    return;
}

bool Verify(ZZX Kv,ZZX s[2], vec_ZZ msg){
    //Hash id
    const ZZX phi = Cyclo();
    ZZ_pX aux = conv<ZZ_pX>((conv<ZZX>(msg) - (s[1]*Kv))%phi);

    s[0] = conv<ZZX>(aux);

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
    return (norm3 < conv<ZZ>(1.36*q0/2*sqrt(2*N0)));
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

        vec_ZZ msg = RandomVector();
        ZZX s[2];
        t1 = clock();
        Sign(s,msg,MSKD);
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
        bool valid = Verify(conv<ZZX>(Kv),s, msg);
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
        cout << "KEMKeyGen: " << t_keygen << "\n";
        cout << "KEMEnc   : " << t_sig << "\n";
        cout << "KEMDec   : " << t_ver << "\n";
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



// void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK)
// {
//     MPKD->h = MPK;
//     ZZXToFFT(MPKD->h_FFT, conv<ZZX>(MPK));
// }



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


unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id, const MSK_Data * const MSKD)
{
    unsigned int i;
    ZZX f,g,t,aux;

    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    
    t = conv<ZZX>(id);
    aux = ((SK_id[0] - t)*f + g*SK_id[1])%phi;

    for(i=0; i<N0; i++)
    {
        aux[i] %= q1;
    }

    if( IsZero(aux) != 0)
    {
        cout << "The signature (s1,s2) doesn't verify the required equality [ (s1 - t)*f + g*s2 = 0 ] !\nActually, (s1 - t)*f + g*s2 = " << aux << endl << endl;
    }
    return IsZero(aux);
}


// void IBE_Encrypt(long C[2][N0], const long m[N0], const long id0[N0], const MPK_Data * const MPKD)
// {

//     unsigned long i;
//     long r[N0], e1[N0], e2[N0];
//     CC_t r_FFT[N0], t_FFT[N0], aux1_FFT[N0], aux2_FFT[N0];

//     for(i=0; i<N0; i++)
//     {
//         e1[i] = (rand()%3) - 1;
//         e2[i] = (rand()%3) - 1;
//         r[i] = (rand()%3) - 1;
//     }

//     MyIntFFT(r_FFT, r);
//     MyIntFFT(t_FFT, id0);

//     for(i=0; i<N0; i++)
//     {
//         aux1_FFT[i] = r_FFT[i]*((MPKD->h_FFT)[i]);
//         aux2_FFT[i] = r_FFT[i]*t_FFT[i];
//     }

//     MyIntReverseFFT(C[0], aux1_FFT);
//     MyIntReverseFFT(C[1], aux2_FFT);

//     for(i=0; i<N0; i++)
//     {
//         C[0][i] = (C[0][i] + e1[i]               + q0/2)%q0 - (q0/2);
//         C[1][i] = (C[1][i] + e2[i] + (q0/2)*m[i] + q0/2)%q0 - (q0/2);
//     } 

// }


// void IBE_Decrypt(long message[N0], const long C[2][N0], const CC_t * const SKid_FFT)
// {
//     unsigned int i;
//     CC_t c0_FFT[N0], aux_FFT[N0];

//     MyIntFFT(c0_FFT, C[0]);

//     for(i=0; i<N0; i++)
//     {
//         aux_FFT[i] = c0_FFT[i]*SKid_FFT[i];
//     }

//     MyIntReverseFFT(message, aux_FFT);

//     for(i=0; i<N0; i++)
//     {
//         message[i] = C[1][i] - message[i];
//         message[i] = ((unsigned long)(message[i] ))%q0;
//         message[i] = (message[i] + (q0>>2) )/(q0>>1);
//         message[i] %= 2;
//     }

// }



//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================


// void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD)
// {
//     clock_t t1, t2;
//     float diff;
//     unsigned int i;
//     vec_ZZ id;
//     ZZX SK_id[2];

//     t1 = clock();

//     cout << "0%" << flush;
//     for(i=0; i<nb_extr; i++)
//     {
//         id = RandomVector();

//         IBE_Extract(SK_id, id, MSKD);
//         if((i+1)%(nb_extr/10)==0)
//         {
//             cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
//         }
//     }

//     t2 = clock();
//     diff = ((float)t2 - (float)t1)/1000000.0F;
//     cout << "\n\nIt took " << diff << " seconds to extract " << nb_extr << " keys." << endl;
//     cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per key." << endl << endl;
// }


// void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
// {
//     clock_t t1, t2;
//     double diff;
//     unsigned int i,j;
//     vec_ZZ id;
//     ZZX SK_id[2], w;
//     CC_t SKid_FFT[N0];
//     long int message[N0], decrypted[N0];
//     long int identity[N0], Ciphertext[2][N0];


//     id = RandomVector();
//     IBE_Extract(SK_id, id, MSKD);
//     IBE_Verify_Key(SK_id, id, MSKD);
//     ZZXToFFT(SKid_FFT, SK_id[1]);
//     for(i=0; i<N0; i++)
//     {
//         identity[i] = conv<long int>(id[i]);
//     }

//     t1 = clock();

//     cout << "0%" << flush;
//     for(i=0; i<nb_cryp; i++)
//     {

//         for(j=0; j<N0; j++)
//         {
//             message[j] = (rand()%2);
//         }

//         IBE_Encrypt(Ciphertext, message, identity, MPKD);
//         IBE_Decrypt(decrypted, Ciphertext, SKid_FFT);

//         if((i+1)%(nb_cryp/10)==0)
//         {
//             cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
//         }
//     }

//     t2 = clock();
//     diff = ((double)t2 - (double)t1)/1000000.0l;
//     cout << "\n\nIt took " << diff << " seconds to do " << nb_cryp << " encryptions and decryptions." << endl;
//     cout << "That's " << (diff/nb_cryp)*1000 << " milliseconds per encryption+decryption." << endl;
//     cout << "That's " << (diff/nb_cryp)*1000*1024/N0 << " milliseconds per encryption+decryption per Kilobit." << endl << endl;
// }


// void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD)
// {
//     unsigned int i, rep;
//     vec_ZZ id;
//     ZZX SK_id[2];

//     rep = 0;

//     cout << "0%" << flush;
//     for(i=0; i<nb_extr; i++)
//     {
//         id = RandomVector();

//         IBE_Extract(SK_id, id, MSKD);
//         rep += IBE_Verify_Key(SK_id, id, MSKD);
//         if((i+1)%(nb_extr/10)==0)
//         {
//             cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
//         }
//     }

//     cout << endl;
//     if(rep == 0)
//     {    cout << endl << nb_extr << " extractions successfully performed!" << endl << endl;    }
//     else
//     {    cout << endl << rep << " out of " << nb_extr << " extractions failed miserabily!" << endl << endl;    }
// }



// void DigSig2(const unsigned int nb_extr, MSK_Data * MSKD)
// {
//     unsigned int i, rep;
//     vec_ZZ id;
//     ZZX SK_id[2];

//     rep = 0;


//     //Gen Key
//     cout <<"\n--1. Keygen--\n Done before this actual point.";
//         //Build MSKD (done before on IBE.cc)

//     //Sig
//     cout <<"\n--2. Signature--";
//         //Using id as our random vector
    
//     id = RandomVector();

//     IBE_Extract(SK_id, id, MSKD);
//     rep += IBE_Verify_Key(SK_id, id, MSKD);

    
//     ZZX f = MSKD -> PrK[0];
//     ZZX g = MSKD -> PrK[1];
//     ZZX t = conv<ZZX>(id);



//     cout << "\nt: ";
//     for (i=0; i < N0; i++){
//         cout << t[i] << ",";
//     }

//     cout << "\nf: ";
//     for (i=0; i < N0; i++){
//         cout << f[i] << ",";
//     }

//     cout << "\ng: ";
//     for (i=0; i < N0; i++){
//         cout << g[i] << ",";
//     }

//     cout << "\nS1: ";
//     for (i=0; i < N0; i++){
//         cout << SK_id[0][i] << ",";
//     }

//     cout << "\nS2: ";
//     for (i=0; i < N0; i++){
//         cout << SK_id[1][i] << ",";
//     }

//     //Verify

//     cout << "\n--Verify---\n";
//     //T has been conv'ed above.
//     ZZX h = conv<ZZX>(Quotient(f,g));
  
//     cout << "\nh: ";
//     for (i=0; i < N0; i++){
//         cout << h[i] << ",";
//     }



//     ZZX aux = (t - (SK_id[1]*h))%phi;

//     cout << "\nRecovered S1: ";
//     for (i=0; i < N0; i++){
//         cout << aux[i]%q0 << ",";
//     }

//     //aux = ((SK_id[0] - t)*f + g*SK_id[1])%phi;

//     cout <<"\n *Mosca (f and g) & (s1 and s2) are used oppositly than our paper. Also notice the ordering is in accending as per usual.";

//     cout << endl;
// }


// void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
// {
//     unsigned int i, j, rep;
//     vec_ZZ id;
//     ZZX SK_id[2], m, w;
//     CC_t SKid_FFT[N0];
//     long int id0[N0], Ciphertext[2][N0];
//     long int message[N0], decrypted[N0];


//     id = RandomVector();
//     IBE_Extract(SK_id, id, MSKD);
//     IBE_Verify_Key(SK_id, id, MSKD);
//     ZZXToFFT(SKid_FFT, SK_id[1]);

//     rep = 0;

//     for(i=0; i<N0; i++)
//     {
//         id0[i] = conv<long int>(id[i]);
//     }

//     cout << "0%" << flush;
//     for(i=0; i<nb_cryp; i++)
//     {

//         for(j=0; j<N0; j++)
//         {
//             message[j] = (rand()%2);
//         }

//         IBE_Encrypt(Ciphertext, message, id0, MPKD);
//         IBE_Decrypt(decrypted, Ciphertext, SKid_FFT);
        
//         for(j=0; j<N0; j++)
//         {
//             if(message[j] != decrypted[j])
//             {
//                 cout << "ERROR : Dec(Enc(m)) != m " << endl;
//                 rep++;
//                 break;
//             }
//         }

//         if((i+1)%(nb_cryp/10)==0)
//         {
//             cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
//         }
//     }

//     cout << endl;
//     if(rep == 0)
//     {    cout << endl << nb_cryp << " encryptions+decryptions successfully performed!" << endl << endl;    }
//     else
//     {    cout << endl << rep << " out of " << nb_cryp << " encryptions+decryptions failed miserabily!" << endl << endl;    }
// }
