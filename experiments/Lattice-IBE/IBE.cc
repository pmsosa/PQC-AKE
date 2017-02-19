/*

Copyright or © or Copr. Thomas Prest.

Thomas.Prest@ENS.fr

This software is a computer program which purpose is to provide to the 
research community a proof-of-concept implementation of the identity-based
encryption scheme over NTRU lattices, described in the paper
"Efficient Identity-Based Encryption over NTRU Lattices", of
Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at
homepages.cwi.nl/~ducas/ , www.di.ens.fr/~lyubash/
and www.di.ens.fr/~prest/ .

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/



#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>


#include "params.h"
#include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "Scheme.h"


using namespace std;
using namespace NTL;


//////////////// DIGITAL SIGNATURE ////////////////

void InitDS(){

    srand(rdtsc());                 //Seed Random

    ZZ_p::init(q1);                 //q1 = ZZ(q0)

    const ZZX phi = Cyclo();        //Phi - Used all over the place
    ZZ_pX phiq = conv<ZZ_pX>(phi);  
    ZZ_pXModulus PHI(phiq);
}

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
    IBE_Extract(s, msg, MSKD);
    return;
}

bool Verify(ZZX Kv,ZZX s[2], vec_ZZ msg, MSK_Data* MSKD){
    //Hash id
    const ZZX phi = Cyclo();
    ZZ_pX aux = conv<ZZ_pX>((conv<ZZX>(msg) - (s[1]*Kv))%phi);

    s[0] = conv<ZZX>(aux);
    return true;
}

void run_DS_example(){

    /// RUNNING THE KEYGEN ALGORITHM
        cout << "\n---Keygen---\n";
        ZZX Ks[2];
        ZZ_pX Kv;
        MSK_Data * MSKD = new MSK_Data;
        SigKeyGen(Ks,Kv,MSKD);
        //ZZX Kv2 = conv<ZZX>(Kv);

        //PRINTIN DEBUG INFO 
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

    /// RUNNING THE SIGN ALGORITHM
        cout << "\n---Sign---\n";

        vec_ZZ msg = RandomVector();
        ZZX s[2];
        Sign(s,msg,MSKD);

        //PRINTIN DEBUG INFO
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

    /// RUNNING THE VERIFY ALGORITHM
        cout << "\n--Verify---\n";

        //Notice
        //1. Turn Kv into ZZX(Kv)
        //2. We will pretend that s[0] is empty since we didn't get that info
        bool valid = Verify(conv<ZZX>(Kv),s, msg, MSKD);

        //PRINTIN DEBUG INFO
        cout << "     derived s1: ";
        for (int i=0; i < N0; i++){
            cout << s[0][i] << ",";
        }
        cout <<"\n";

        cout << "          Valid: " << valid <<"(not yet implemented) \n" ;

    return;
}


/////////// Key Encapsulation Mechanism ///////////
void KEMKeyGen(ZZX& Kd, ZZX& Ke){
    
    ZZ small_norm = conv<ZZ>(30);

    bool valid = false;
    while (!valid){

        try{
            //Build f & g        
            ZZX f = RandomPolyFixedSqNorm(small_norm,N0-1);
            ZZX g = RandomPolyFixedSqNorm(small_norm,N0-1);
            Kd = g;
            

            //Check that f is invertible in both Zq[x]/<x^n+1> and Z2[x]/<x^n+1>
            ZZX f_inv = conv<ZZX>(Inverse2(conv<ZZX>(f),q0));
            ZZX f_inv2 = conv<ZZX>(Inverse2(conv<ZZX>(f),2));

            //Check that g is invertible in both Zq[x]/<x^n+1> and Z2[x]/<x^n+1>
            ZZX g_inv = conv<ZZX>(Inverse2(conv<ZZX>(g),q0));
            ZZX g_inv2 = conv<ZZX>(Inverse2(conv<ZZX>(g),2));


            //Debuggin' Information
            cout << "f: "<< f << " ("<<deg(f)<<")"<<"\n";
            cout << "g: "<< g << " ("<<deg(g)<<")"<<"\n";

            cout << "f^-1: "<< f_inv << " ("<<deg(f_inv)<<")"<<"\n";
            cout << "f^-1 (2): "<< f_inv2 << " ("<<deg(f_inv2)<<")"<<"\n";

            cout << "g^-1: "<< g_inv << " ("<<deg(g_inv)<<")"<<"\n";
            cout << "g^-1 (2): "<< g_inv2 << " ("<<deg(g_inv2)<<")"<<"\n";
            valid = true;

            Ke = (f*g_inv);
            
        }
        catch (int e){
            cout << "Caught!\n";
            //Do nothing!
        }
    }



    return;
}


void Encapsulate(){
    return;
}

void Decapsulate(){
    return;
}

int main(){

    srand(rdtsc());

    //Test the DS
    //run_DS_example();

    //Test the KEM
    ZZX Kd,Ke;
    KEMKeyGen(Kd,Ke);
    cout << Kd << " | " << Ke << "\n";


    return 0;
}









//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================


// void DigSig()
// {

//     cout << "\n====Digital Signature====\n";

//     /////////////////////////////////////////////////////////1. KEYGEN
//     cout <<"\n--1. Keygen--\n";

//     ZZ_p::init(q1);
//     zz_p::init(q0);
    
//     ZZX MSK[4];
//     ZZ_pX phiq, MPK;
//     MSK_Data * MSKD = new MSK_Data;
//     //MPK_Data * MPKD = new MPK_Data;
//     const ZZX phi = Cyclo();
//     srand(rdtsc()); // initialisation of rand

//     phiq = conv<ZZ_pX>(phi);
//     ZZ_pXModulus PHI(phiq);
//     Keygen(MPK, MSK);
//     CompleteMSK(MSKD, MSK);
//     //CompleteMPK(MPKD, MPK);


//     /////////////////////////////////////////////////////////2. Sign Message
//     cout <<"\n--2. Signature--";
//         //Using id as our random vector

//     unsigned int i;
//     vec_ZZ id;
//     ZZX SK_id[2];

//     id = RandomVector();



//     IBE_Extract(SK_id, id, MSKD);
//     //IBE_Verify_Key(SK_id, id, MSKD);

    
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
//     cout << "  (";
//     for (i=0; i < N0; i++){
//         cout << SK_id[0][i]%q0 << ",";
//     }
//     cout << " )";

//     cout << "\nS2: ";
//     for (i=0; i < N0; i++){
//         cout << SK_id[1][i] << ",";
//     }

//     /////////////////////////////////////////////////////////3. Verify Message

//     cout << "\n--3. Verify---\n";
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

// int main()
// {


//     DigSig();

//     return 0;
//     cout << "\n=======================================================================\n";
//     cout << "This program is a proof-of concept for efficient IBE over lattices.\n";
//     cout << "It generates a NTRU lattice of dimension 2N and associated modulus q,\n";
//     cout << "and perform benches and tests, for user key extraction and encryption/decryption.";
//     cout << "\n=======================================================================\n\n";

//     ZZX MSK[4];
//     ZZ_pX phiq, MPK;
//     unsigned int i;
//     float diff;
//     MSK_Data * MSKD = new MSK_Data;
//     MPK_Data * MPKD = new MPK_Data;
//     clock_t t1, t2;
//     const ZZX phi = Cyclo();

//     srand(rdtsc()); // initialisation of rand

//     cout << "N = " << N0 << endl;
//     cout << "q = " << q0 << endl;

//     ZZ_p::init(q1);
//     //zz_p::init(q0);

//     phiq = conv<ZZ_pX>(phi);
//     ZZ_pXModulus PHI(phiq);


//     cout << "\n===================================================================\n KEY GENERATION";
//     cout << "\n===================================================================\n";
//     t1 = clock();
//     for(i=0; i<1; i++)
//     {
//         Keygen(MPK, MSK);
//     }

//     CompleteMSK(MSKD, MSK);
//     CompleteMPK(MPKD, MPK);

//     t2 = clock();
//     diff = ((float)t2 - (float)t1)/1000000.0F;
//     cout << "It took " << diff << " seconds to generate the Master Secret Key" << endl;



//     // //==============================================================================
//     // //Key extraction bench and encryption/decryption bench
//     // //==============================================================================
//     // const unsigned int nb_extrb = 100;
//     // const unsigned int nb_crypb = 1000;

//     // cout << "\n===================================================================\n RUNNING EXTRACTION BENCH FOR ";
//     // cout << nb_extrb << " DIFFERENT IDENTITIES\n===================================================================\n";
//     // Extract_Bench(nb_extrb, MSKD);

//     // cout << "\n===================================================================\n RUNNING ENCRYPTION BENCH FOR ";
//     // cout << nb_crypb << " DIFFERENT MESSAGES\n===================================================================\n";
//     // Encrypt_Bench(nb_crypb, MPKD, MSKD);



//     ///==============================================================================
//     //Key extraction test and encryption/decryption test
//     //==============================================================================
//     const unsigned int nb_extrt = 100;
//     const unsigned int nb_crypt = 100;

//     cout << "\n===================================================================\n CHECKING EXTRACTION VALIDITY FOR ";
//     cout << nb_extrt << " DIFFERENT IDENTITIES\n===================================================================\n";
//     Extract_Test(nb_extrt, MSKD);

//     cout << "\n=============DIGSIG============\n";
//     DigSig2(nb_extrt, MSKD);

//     // cout << "\n===================================================================\n CHECKING ENCRYPTION VALIDITY FOR ";
//     // cout << nb_extrt << " DIFFERENT MESSAGES\n===================================================================\n";
//     // Encrypt_Test(nb_crypt, MPKD, MSKD);

//     free(MSKD);
//     free(MPKD);
//     return 0;
// }
