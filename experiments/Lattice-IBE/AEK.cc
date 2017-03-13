////
// AEK Implementation based on Del Pino, Lyubasgevsky, Pointcheval's work.
// Some of the code has been repurposed from Thomas.Prest@ENS.fr
// Pedro M. Sosa
////

//The Following is the Copyright from Thomas Prest's code:
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
#include "DigitalSignature.h"
#include "KEM.h"

#include <openssl/sha.h>



using namespace std;
using namespace NTL;

//const ZZX phi = Cyclo();



void Hash2(vec_ZZ& Auth, vec_ZZ Ke,vec_ZZ c,vec_ZZ k){
    SHA256_CTX ctx;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Init(&ctx);

    //vec_ZZ msg = RandomVector();
    Auth = Ke;
    for(int i = 0; i < N0; i++) {
        char f = (conv<int>(Ke[i]+c[i]+k[i])%255);
        SHA256_Update(&ctx, &f, 1);
    }
    SHA256_Final(digest, &ctx);

    for(int i=0; i < SHA256_DIGEST_LENGTH and i < N0; i++){
        //cout << int(digest[i]) << " ";
        Auth[i] = conv<ZZ>(digest[i])%q0;
    }

}


void Hash1(vec_ZZ& sk, vec_ZZ Ke,vec_ZZ c, vec_ZZ auth, vec_ZZ k){
    SHA256_CTX ctx;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Init(&ctx);
    SHA256_Update(&ctx, "K", 1);

    //vec_ZZ msg = RandomVector();
    sk = Ke;
    for(int i = 0; i < N0; i++) {
        char f = (conv<int>(Ke[i]+c[i]+k[i]+auth[i])%255);
        SHA256_Update(&ctx, &f, 1);
    }
    SHA256_Final(digest, &ctx);

    for(int i=0; i < SHA256_DIGEST_LENGTH and i < N0; i++){
        //cout << int(digest[i]) << " ";
        sk[i] = conv<ZZ>(digest[i])%q0;
    }
}


//Sometimes all you need is to take a little break.
void BREAK(int message){
	int a = 0;
	cout << "Break!"<< message <<"\n";
	cin >> a;
}

void AKE_example(){

	//Alice and Bob Generate SigKeyGen (Previous to actual AKE Exchange)
		ZZX 	Ks_a[2], Ks_b[2];
		ZZ_pX 	Kv_a   , Kv_b;

		MSK_Data * MSKD_a = new MSK_Data;													//ALice: (Ks1,Kv1) <- SigKeyGen
		MSK_Data * MSKD_b = new MSK_Data;													//Bob:   (Ks2,Ks2) <- SigKeyGen

		SigKeyGen(Ks_a,Kv_a,MSKD_a);														//ALICE: sk(1) <- BOT
		SigKeyGen(Ks_b,Kv_b,MSKD_b);														//BOB:   sk(2) <- BOT



	//Alice
		vec_ZZ sk_a;
		ZZX Kd_a, Ke_a_temp;
		KEMKeyGen(Kd_a,Ke_a_temp);															//ALICE: (Kd,Ke) <- KEMKeyGen
		vec_ZZ Ke_a = conv<vec_ZZ>(Ke_a_temp);

		cout << "Ke_a: "<< Ke_a <<"\n";

		vec_ZZ r_a;
		ZZX s_a[2];

		Sign(s_a,Ke_a,r_a,MSKD_a); //Ks_a == MSKD_a which contains f,g 						//ALICE: simga1 <- Sig(Ks1,Ke)
								   //s1 => renamed as s2 will be inside s_a[1]

	//Alice sends Sigma1 to Bob 															//AlICE ----> simga1 = <Ke>1 ----> Bob
		// Aka. Bob has access only to:
		// - s_a[1] ]
		// - Ke_a   ]-- Sent from Alice
		// - r_a    ]
		// - Kv_a   ]--- (Public) Obtained before hand.


	//Bob
		vec_ZZ sk_b;
		ZZX c_b,k_b; 
		vec_ZZ Auth_b;
		ZZX s_b[2];	
		vec_ZZ r_b; 
		//cout << Kv_a;
		if ( Verify(conv<ZZX>(Kv_a),s_a,Ke_a,r_a) ){										//BOB: if (Ver(Kv1,simga1) != BOT) Then
			
			
			Encapsulate(Ke_a_temp,c_b,k_b);													//BOB: (c,k) <- Enc(Ke)
			
			Hash2(Auth_b, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_b));						//BOB: Auth  <- H2(sigma1,c,k)

			vec_ZZ c_Auth;
			c_Auth = Auth_b;
			c_Auth.append(conv<vec_ZZ>(c_b));

			Sign(s_b, Auth_b, r_b, MSKD_b);	//Ks_a == MSKD_a which contains f,g 													//BOB: Sigma2 <- Sig(Ks2,(c,k))
											//s1 => renamed as s2 will be inside s_a[1]

			Hash1(sk_b, Ke_a, conv<vec_ZZ>(c_b), Auth_b, conv<vec_ZZ>(k_b));
			cout <<"SK_b = "<<sk_b <<"\n";
		}
		else{ cout << "Bob: Abort!"; return; }


	//Bob sends Sigma2 to Alice
		// Aka. Alice has access to:
		// c_b    ]
		// Auth_b ] -- Sent from Bob
		// Kv_b   ] -- (Public Obtained before hand).

	//Alice
		if ( Verify(conv<ZZX>(Kv_b),s_b,Auth_b,r_b) ){
			ZZX k_a;
			Decapsulate(Kd_a,c_b,k_a);

			vec_ZZ Auth_a;
			Hash2(Auth_a, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_a));
			if (IsZero(Auth_a - Auth_b)){
				Hash1(sk_a, Ke_a, conv<vec_ZZ>(c_b), Auth_b, conv<vec_ZZ>(k_b));
				cout <<"SK_a = "<<sk_a <<"\n";
				//sk_1 <- H1(sigma1,sigma2,k') 
			}
			else{cout << "Alice: Abort! (H2 Hashes Didn't Match!)";}
			

		}
		else{ cout << "Alice: Abort!"; return; }	


	//CHECK IF EVERYTHING WORKED!
		if (IsZero(sk_a - sk_b)){
			cout << "\nsk_a == sk_b; Successful AEK!\n";
		}


}

int main(){

	

    srand(rdtsc());

    if (false){
	    cout <<"\n\n--RUNNING THE DS/KEM Examples--\n\n";
	    //Test the DS
	    run_DS_example();

	    //Test the KEM
	    run_KEM_example();
	}	
	
    if (true){
	    //AKE Example
	    cout <<"\n\n--RUNNING THE AEK Examples--\n\n";
	    AKE_example();
	}
    return 0;
}





