////
// AKE Implementation based on Del Pino, Lyubasgevsky, Pointcheval's work.
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
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "DigitalSignature.h"
#include "KEM.h"
#include "cpucycles.h"

#include <openssl/sha.h>
#include "huffman.h"




using namespace std;
using namespace NTL;

//const ZZX phi = Cyclo();



void Hash2(vec_ZZ& Auth, vec_ZZ Ke,vec_ZZ c,vec_ZZ k){
    SHA256_CTX ctx;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Init(&ctx);

    //vec_ZZ msg = RandomVector();
    Auth.SetLength(32);
    for(int i = 0; i < N0; i++) {
        char f = (conv<int>(Ke[i]+c[i]+k[i])%255);
        SHA256_Update(&ctx, &f, 1);
    }
    SHA256_Final(digest, &ctx);

    for(int i=0; i < 32 and i < N0; i++){
        //cout << int(digest[i]) << " ";
        Auth[i] = conv<ZZ>(digest[i])%q0;
    }

}


void Hash1(vec_ZZ& sk, vec_ZZ Ke,vec_ZZ c_Auth, vec_ZZ k){
    SHA256_CTX ctx;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Init(&ctx);
    SHA256_Update(&ctx, "K", 1);

    //vec_ZZ msg = RandomVector();
    sk.SetLength(32);
    for(int i = 0; i < N0; i++) {
        char f = (conv<int>(Ke[i]+c_Auth[i]+k[i])%255);
        SHA256_Update(&ctx, &f, 1);
    }
    SHA256_Final(digest, &ctx);

    for(int i=0; i < 32 and i < N0; i++){
        //cout << int(digest[i]) << " ";
        sk[i] = conv<ZZ>(digest[i])%q0;
    }
}


void AKE_timed_example(){

	clock_t t1,t2;
	float ta_ds_keygen, ta_kem_keygen, ta_ds_sig, ta_ds_ver, ta_kem_dec, ta_h2, ta_h1;
	float tb_ds_keygen, tb_ds_ver, tb_kem_enc, tb_h2, tb_ds_sig, tb_h1;

	//Alice and Bob Generate SigKeyGen (Previous to actual AKE Exchange)
		ZZX 	Ks_a[2], Ks_b[2];
		ZZ_pX 	Kv_a   , Kv_b;

		MSK_Data * MSKD_a = new MSK_Data;													//ALice: (Ks1,Kv1) <- SigKeyGen
		MSK_Data * MSKD_b = new MSK_Data;													//Bob:   (Ks2,Ks2) <- SigKeyGen

			t1 = clock();
		SigKeyGen(Ks_a,Kv_a,MSKD_a);														//ALICE: sk(1) <- BOT
			t2 = clock();
			ta_ds_keygen = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;
		


			t1 = clock();
		SigKeyGen(Ks_b,Kv_b,MSKD_b);														//BOB:   sk(2) <- BOT
			t2 = clock();
			tb_ds_keygen = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;


	//Alice
		vec_ZZ sk_a;
		ZZX Kd_a, Ke_a_temp, Kd_inv2_a;

			t1 = clock();
		KEMKeyGen(Kd_a,Ke_a_temp, Kd_inv2_a);															//ALICE: (Kd,Ke) <- KEMKeyGen
			t2 = clock();
			ta_kem_keygen = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

		vec_ZZ Ke_a = conv<vec_ZZ>(Ke_a_temp);

		//cout << "Ke_a: "<< Ke_a <<"\n";

		vec_ZZ r_a;
		ZZX s_a[2];

			t1 = clock();
		Sign(s_a,Ke_a,r_a,MSKD_a); //Ks_a == MSKD_a which contains f,g 						//ALICE: simga1 <- Sig(Ks1,Ke)
								   //s1 => renamed as s2 will be inside s_a[1]

			t2 = clock();
			ta_ds_sig = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			//cout << "S:" << deg(s_a[1]) << Ke_a<< "\n";

	//Alice sends Sigma1 to Bob 															//AlICE ----> simga1 = <Ke>1 ----> Bob
		// Aka. Bob has access only to:
		// - s_a[1] ] 
		// - Ke_a   ]  -- Sent from Alice
		// - r_a    ] 
		// - Kv_a   ]--- (Public) Obtained before hand.
		cout << "\nTransmitting (A->b)\n";
		cout << " +    Ke_a : " << Ke_a.length()*16 << "\n";
		cout << " +      s1 : " << conv<vec_ZZ>(s_a[1]).length()*8 << "\n";
		cout << " +     r_b : " << r_a.length()*8 << "\n";
		cout << "Total bits : " << r_a.length()*8+Ke_a.length()*16+conv<vec_ZZ>(s_a[1]).length()*8 << "\n";


	//Bob
		vec_ZZ sk_b;
		ZZX c_b,k_b; 
		vec_ZZ Auth_b;
		ZZX s_b[2];	
		vec_ZZ r_b; 
		vec_ZZ c_Auth;
		//cout << Kv_a;
			t1 = clock();
		if ( Verify(conv<ZZX>(Kv_a),s_a,Ke_a,r_a) ){										//BOB: if (Ver(Kv1,simga1) != BOT) Then
				t2 = clock();
				tb_ds_ver = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;
			
				t1 = clock();
			Encapsulate(Ke_a_temp,c_b,k_b);													//BOB: (c,k) <- Enc(Ke)
				t2 = clock();
				tb_kem_enc = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

				t1 = clock();
			Hash2(Auth_b, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_b));						//BOB: Auth  <- H2(sigma1,c,k)
				t2 = clock();
				tb_h2 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			
			c_Auth = Auth_b;
			c_Auth.append(conv<vec_ZZ>(c_b));

			

				t1 = clock();
			Sign(s_b, c_Auth, r_b, MSKD_b);	//Ks_a == MSKD_a which contains f,g 													//BOB: Sigma2 <- Sig(Ks2,(c,k))
											//s1 => renamed as s2 will be inside s_a[1]



				t2 = clock();
				tb_ds_sig = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

				t1 = clock();
			Hash1(sk_b, Ke_a, c_Auth, conv<vec_ZZ>(k_b));
				t2 = clock();
				tb_h1 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			
			//cout <<"SK_b = "<<sk_b <<"\n";
		}
		else{ cout << "Bob: Abort!"; return; }


	//Bob sends Sigma2 to Alice
		// Aka. Alice has access to:
		// c_Auth ] 
		// s1_b   ] -- Sent from Bob
		// r_b    ]
		// Kv_b   ] -- (Public Obtained before hand).
		cout << "\nTransmitting (B->A)\n";
		cout << " +  c_auth : " << c_Auth.length()*16 << "\n";
		cout << " +      s1 : " << conv<vec_ZZ>(s_b[1]).length()*8 << "\n";
		cout << " +     r_b : " << r_b.length()*8 << "\n";
		cout << "Total bits : " <<conv<vec_ZZ>(s_b[1]).length()*8+c_Auth.length()*16+r_b.length()*8 << "\n";


	//Alice
			t1 = clock();
		if ( Verify(conv<ZZX>(Kv_b),s_b,c_Auth,r_b) ){
				t2 = clock();
				ta_ds_ver = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			ZZX k_a;
			
				t1 = clock();
			Decapsulate(Kd_a,c_b,k_a, Kd_inv2_a);
				t2 = clock();
				ta_kem_dec = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			vec_ZZ Auth_a;
			
				t1 = clock();
			Hash2(Auth_a, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_a));
				t2 = clock();
				ta_h2 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			if (IsZero(Auth_a - Auth_b)){

					t1 = clock();
				Hash1(sk_a, Ke_a, c_Auth, conv<vec_ZZ>(k_b));
					t2 = clock();
					ta_h1 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

				//cout <<"SK_a = "<<sk_a <<"\n";

			}
			else{cout << "Alice: Abort! (H2 Hashes Didn't Match!)"; return; }
			

		}
		else{ cout << "Alice: Abort!"; return; }	


	//CHECK IF EVERYTHING WORKED!
		if (IsZero(sk_a - sk_b)){
			cout << "\nsk_a == sk_b; Successful AKE!\n";
		}


		cout << "N=" << N0 <<", q=" << q0 << "\n";
		cout << "KEM_Norm=" << KEM_NORM << ", DS_Sigma=" << DS_SIGMA <<"\n\n";
		cout << "\n\n========Time (ms)=======\n";
		cout << "\n==ALICE==\n";
		cout << "- DS_keygen  :"<< ta_ds_keygen		<< "\n\n";
		cout << "- KEM_keygen :"<< ta_kem_keygen	<< "\n";
		cout << "- DS_sign    :"<< ta_ds_sig		<< "\n";
		cout << "- DS_verify  :"<< ta_ds_ver		<< "\n";
		cout << "- KEM_dec    :"<< ta_kem_dec		<< "\n";
		cout << "- H2         :"<< ta_h2			<< "\n";
		cout << "- H1         :"<< ta_h1			<< "\n";
		cout << " TOTAL (ALL) :"<< ta_ds_keygen + ta_kem_keygen + ta_ds_sig + ta_ds_ver + ta_kem_dec + ta_h2 + ta_h1 << "\n";
		cout << " TOTAL (AKE) :"<< ta_kem_keygen + ta_ds_sig + ta_ds_ver + ta_kem_dec + ta_h2 + ta_h1 << "\n";

		cout << "\n==BOB==\n";
		cout << "- DS_keygen  :"<< tb_ds_keygen		<< "\n\n";
		cout << "- DS_verify  :"<< tb_ds_ver		<< "\n";
		cout << "- KEM_enc    :"<< tb_kem_enc		<< "\n";
		cout << "- H2         :"<< tb_h2			<< "\n";
		cout << "- DS_sign    :"<< tb_ds_sig		<< "\n";
		cout << "- H1         :"<< tb_h1			<< "\n";
		cout << " TOTAL (ALL) :"<< tb_ds_keygen + tb_ds_ver + tb_kem_enc + tb_h2 + tb_ds_sig + tb_h1 << "\n";
		cout << " TOTAL (AKE) :"<< tb_ds_ver + tb_kem_enc + tb_h2 + tb_ds_sig + tb_h1 << "\n";

}


void AKE_timed_exampleMR(){

	clock_t t1,t2;
	float ta_ds_keygen, ta_kem_keygen, ta_ds_sig, ta_ds_ver, ta_kem_dec, ta_h2, ta_h1;
	float tb_ds_keygen, tb_ds_ver, tb_kem_enc, tb_h2, tb_ds_sig, tb_h1;

	//Alice and Bob Generate SigKeyGen (Previous to actual AKE Exchange)
		ZZX 	Ks_a[2], Ks_b[2];
		ZZ_pX 	Kv_a   , Kv_b;

		MSK_Data * MSKD_a = new MSK_Data;													//ALice: (Ks1,Kv1) <- SigKeyGen
		MSK_Data * MSKD_b = new MSK_Data;													//Bob:   (Ks2,Ks2) <- SigKeyGen

			t1 = clock();
		SigKeyGen(Ks_a,Kv_a,MSKD_a);														//ALICE: sk(1) <- BOT
			t2 = clock();
			ta_ds_keygen = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;
		


			t1 = clock();
		SigKeyGen(Ks_b,Kv_b,MSKD_b);														//BOB:   sk(2) <- BOT
			t2 = clock();
			tb_ds_keygen = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;


	//Alice
		vec_ZZ sk_a;
		ZZX Kd_a, Ke_a_temp, Kd_inv2_a;

			t1 = clock();
		KEMKeyGen(Kd_a,Ke_a_temp, Kd_inv2_a);															//ALICE: (Kd,Ke) <- KEMKeyGen
			t2 = clock();
			ta_kem_keygen = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

		vec_ZZ Ke_a = conv<vec_ZZ>(Ke_a_temp);

		//cout << "Ke_a: "<< Ke_a <<"\n";

		vec_ZZ r_a;
		ZZX s_a[2];
		vec_ZZ m2_a;

			t1 = clock();
		//TODO
		Sign2(s_a,m2_a,Ke_a,MSKD_a); //Ks_a == MSKD_a which contains f,g 						//ALICE: simga1 <- Sig(Ks1,Ke)
								   //s1 => renamed as s2 will be inside s_a[1]

			t2 = clock();
			ta_ds_sig = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			//cout << "S:" << deg(s_a[1]) << Ke_a<< "\n";

	//Alice sends Sigma1 to Bob 															//AlICE ----> simga1 = <Ke>1 ----> Bob
		// Aka. Bob has access only to:
		// - s_a[1] ] 
		// - Ke_a   ]  -- Sent from Alice
		// - r_a    ] 
		// - Kv_a   ]--- (Public) Obtained before hand.
		cout << "\nTransmitting (A->b)\n";
		cout << " +    Ke_a : " << Ke_a.length()*16 << "\n";
		cout << " +      s1 : " << conv<vec_ZZ>(s_a[1]).length()*8 << "\n";
		cout << " +     r_b : " << r_a.length()*8 << "\n";
		cout << "Total bits : " << r_a.length()*8+Ke_a.length()*16+conv<vec_ZZ>(s_a[1]).length()*8 << "\n";


	//Bob
		vec_ZZ sk_b;
		ZZX c_b,k_b; 
		vec_ZZ Auth_b;
		ZZX s_b[2];	
		vec_ZZ r_b; 
		vec_ZZ c_Auth;
		vec_ZZ m1_a;
		vec_ZZ m2_b;
		//cout << Kv_a;
			t1 = clock();
		//TODO
		if ( Verify2(conv<ZZX>(Kv_a),s_a,m2_a,m1_a) ){										//BOB: if (Ver(Kv1,simga1) != BOT) Then
				t2 = clock();
				tb_ds_ver = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;
			Ke_a = m1_a;
			Ke_a.append(m2_a);
			Ke_a_temp = conv<ZZX>(Ke_a);

				t1 = clock();
			Encapsulate(Ke_a_temp,c_b,k_b);													//BOB: (c,k) <- Enc(Ke)
				t2 = clock();
				tb_kem_enc = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

				t1 = clock();
			Hash2(Auth_b, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_b));						//BOB: Auth  <- H2(sigma1,c,k)
				t2 = clock();
				tb_h2 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			
			c_Auth = Auth_b;
			c_Auth.append(conv<vec_ZZ>(c_b));

			
			
				t1 = clock();
			//TODO
			Sign2(s_b, m2_b,c_Auth, MSKD_b);	//Ks_a == MSKD_a which contains f,g 													//BOB: Sigma2 <- Sig(Ks2,(c,k))
											//s1 => renamed as s2 will be inside s_a[1]



				t2 = clock();
				tb_ds_sig = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

				t1 = clock();
			Hash1(sk_b, Ke_a, c_Auth, conv<vec_ZZ>(k_b));
				t2 = clock();
				tb_h1 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			
			//cout <<"SK_b = "<<sk_b <<"\n";
		}
		else{ cout << "Bob: Abort!"; return; }


	//Bob sends Sigma2 to Alice
		// Aka. Alice has access to:
		// c_Auth ] 
		// s1_b   ] -- Sent from Bob
		// r_b    ]
		// Kv_b   ] -- (Public Obtained before hand).
		cout << "\nTransmitting (B->A)\n";
		cout << " +  c_auth : " << c_Auth.length()*16 << "\n";
		cout << " +      s1 : " << conv<vec_ZZ>(s_b[1])*8 << "\n";
		cout << " +     r_b : " << r_b.length()*8 << "\n";
		cout << "Total bits : " <<conv<vec_ZZ>(s_b[1]).length()*8+c_Auth.length()*16+r_b.length()*8 << "\n";

		cout << "s1: "<<s_b[1] <<"\n";

	//Alice
			t1 = clock();
		//TODO

		vec_ZZ m1_b;
		if ( Verify2(conv<ZZX>(Kv_b),s_b,m2_b,m1_b) ){
				t2 = clock();
				ta_ds_ver = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;
			c_Auth = m1_b;
			c_Auth.append(m2_b);

			ZZX k_a;
			
				t1 = clock();
			Decapsulate(Kd_a,c_b,k_a, Kd_inv2_a);
				t2 = clock();
				ta_kem_dec = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			vec_ZZ Auth_a;
			
				t1 = clock();
			Hash2(Auth_a, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_a));
				t2 = clock();
				ta_h2 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

			if (IsZero(Auth_a - Auth_b)){

					t1 = clock();
				Hash1(sk_a, Ke_a, c_Auth, conv<vec_ZZ>(k_b));
					t2 = clock();
					ta_h1 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC * 1000;

				//cout <<"SK_a = "<<sk_a <<"\n";

			}
			else{cout << "Alice: Abort! (H2 Hashes Didn't Match!)"; return; }
			

		}
		else{ cout << "Alice: Abort!"; return; }	


	//CHECK IF EVERYTHING WORKED!
		if (IsZero(sk_a - sk_b)){
			cout << "\nsk_a == sk_b; Successful AKE!\n";
		}


		cout << "N=" << N0 <<", q=" << q0 << "\n";
		cout << "KEM_Norm=" << KEM_NORM << ", DS_Sigma=" << DS_SIGMA <<"\n\n";
		cout << "\n\n========Time (ms)=======\n";
		cout << "\n==ALICE==\n";
		cout << "- DS_keygen  :"<< ta_ds_keygen		<< "\n\n";
		cout << "- KEM_keygen :"<< ta_kem_keygen	<< "\n";
		cout << "- DS_sign    :"<< ta_ds_sig		<< "\n";
		cout << "- DS_verify  :"<< ta_ds_ver		<< "\n";
		cout << "- KEM_dec    :"<< ta_kem_dec		<< "\n";
		cout << "- H2         :"<< ta_h2			<< "\n";
		cout << "- H1         :"<< ta_h1			<< "\n";
		cout << " TOTAL (ALL) :"<< ta_ds_keygen + ta_kem_keygen + ta_ds_sig + ta_ds_ver + ta_kem_dec + ta_h2 + ta_h1 << "\n";
		cout << " TOTAL (AKE) :"<< ta_kem_keygen + ta_ds_sig + ta_ds_ver + ta_kem_dec + ta_h2 + ta_h1 << "\n";

		cout << "\n==BOB==\n";
		cout << "- DS_keygen  :"<< tb_ds_keygen		<< "\n\n";
		cout << "- DS_verify  :"<< tb_ds_ver		<< "\n";
		cout << "- KEM_enc    :"<< tb_kem_enc		<< "\n";
		cout << "- H2         :"<< tb_h2			<< "\n";
		cout << "- DS_sign    :"<< tb_ds_sig		<< "\n";
		cout << "- H1         :"<< tb_h1			<< "\n";
		cout << " TOTAL (ALL) :"<< tb_ds_keygen + tb_ds_ver + tb_kem_enc + tb_h2 + tb_ds_sig + tb_h1 << "\n";
		cout << " TOTAL (AKE) :"<< tb_ds_ver + tb_kem_enc + tb_h2 + tb_ds_sig + tb_h1 << "\n";

}


// void AKE_clocked_example(){

// 	unsigned long long t1,t2;
// 	long long ta_ds_keygen, ta_kem_keygen, ta_ds_sig, ta_ds_ver, ta_kem_dec, ta_h2, ta_h1;
// 	long long tb_ds_keygen, tb_ds_ver, tb_kem_enc, tb_h2, tb_ds_sig, tb_h1;

// 	//Alice and Bob Generate SigKeyGen (Previous to actual AKE Exchange)
// 		ZZX 	Ks_a[2], Ks_b[2];
// 		ZZ_pX 	Kv_a   , Kv_b;

// 		MSK_Data * MSKD_a = new MSK_Data;													//ALice: (Ks1,Kv1) <- SigKeyGen
// 		MSK_Data * MSKD_b = new MSK_Data;													//Bob:   (Ks2,Ks2) <- SigKeyGen

// 			t1 = cpucycles();
// 		SigKeyGen(Ks_a,Kv_a,MSKD_a);														//ALICE: sk(1) <- BOT
// 			t2 = cpucycles();
// 			ta_ds_keygen = t2-t1;
		
// 			t1 = cpucycles();
// 		SigKeyGen(Ks_b,Kv_b,MSKD_b);														//BOB:   sk(2) <- BOT
// 			t2 = cpucycles();
// 			tb_ds_keygen = t2-t1;


// 	//Alice
// 		vec_ZZ sk_a;
// 		ZZX Kd_a, Ke_a_temp, Kd_inv2_a;

// 			t1 = cpucycles();
// 		KEMKeyGen(Kd_a,Ke_a_temp, Kd_inv2_a);															//ALICE: (Kd,Ke) <- KEMKeyGen
// 			t2 = cpucycles();
// 			ta_kem_keygen = t2-t1;

// 		vec_ZZ Ke_a = conv<vec_ZZ>(Ke_a_temp);

// 		//cout << "Ke_a: "<< Ke_a <<"\n";

// 		vec_ZZ r_a;
// 		ZZX s_a[2];

// 			t1 = cpucycles();
// 		Sign(s_a,Ke_a,r_a,MSKD_a); //Ks_a == MSKD_a which contains f,g 						//ALICE: simga1 <- Sig(Ks1,Ke)
// 								   //s1 => renamed as s2 will be inside s_a[1]
// 			t2 = cpucycles();
// 			ta_ds_sig = t2-t1;

// 	//Alice sends Sigma1 to Bob 															//AlICE ----> simga1 = <Ke>1 ----> Bob
// 		// Aka. Bob has access only to:
// 		// - s_a[1] ]
// 		// - Ke_a   ]-- Sent from Alice
// 		// - r_a    ]
// 		// - Kv_a   ]--- (Public) Obtained before hand.

// 	//Bob
// 		vec_ZZ sk_b;
// 		ZZX c_b,k_b; 
// 		vec_ZZ Auth_b;
// 		ZZX s_b[2];	
// 		vec_ZZ r_b; 
// 		//cout << Kv_a;
// 			t1 = cpucycles();
// 		if ( Verify(conv<ZZX>(Kv_a),s_a,Ke_a,r_a) ){										//BOB: if (Ver(Kv1,simga1) != BOT) Then
// 				t2 = cpucycles();
// 				tb_ds_ver = t2-t1;
			
// 				t1 = cpucycles();
// 			Encapsulate(Ke_a_temp,c_b,k_b);													//BOB: (c,k) <- Enc(Ke)
// 				t2 = cpucycles();
// 				tb_kem_enc = t2-t1;

// 				t1 = cpucycles();
// 			Hash2(Auth_b, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_b));						//BOB: Auth  <- H2(sigma1,c,k)
// 				t2 = cpucycles();
// 				tb_h2 = t2-t1;

// 			vec_ZZ c_Auth;
// 			c_Auth = Auth_b;
// 			c_Auth.append(conv<vec_ZZ>(c_b));

// 			cout << "C:" << c_Auth << c_Auth.length() << "\n";

// 				t1 = cpucycles();
// 			Sign(s_b, Auth_b, r_b, MSKD_b);	//Ks_a == MSKD_a which contains f,g 													//BOB: Sigma2 <- Sig(Ks2,(c,k))
// 											//s1 => renamed as s2 will be inside s_a[1]
// 				t2 = cpucycles();
// 				tb_ds_sig = t2-t1;

// 				t1 = cpucycles();
// 			Hash1(sk_b, Ke_a, conv<vec_ZZ>(c_b), Auth_b, conv<vec_ZZ>(k_b));
// 				t2 = cpucycles();
// 				tb_h1 = t2-t1;
			
// 			//cout <<"SK_b = "<<sk_b <<"\n";
// 		}
// 		else{ cout << "Bob: Abort!"; return; }


// 	//Bob sends Sigma2 to Alice
// 		// Aka. Alice has access to:
// 		// c_b    ]
// 		// Auth_b ] -- Sent from Bob
// 		// Kv_b   ] -- (Public Obtained before hand).


// 	//Alice
// 			t1 = cpucycles();
// 		if ( Verify(conv<ZZX>(Kv_b),s_b,Auth_b,r_b) ){
// 				t2 = cpucycles();
// 				ta_ds_ver = t2-t1;

// 			ZZX k_a;
			
// 				t1 = cpucycles();
// 			Decapsulate(Kd_a,c_b,k_a, Kd_inv2_a);
// 				t2 = cpucycles();
// 				ta_kem_dec = t2-t1;

// 			vec_ZZ Auth_a;
			
// 				t1 = cpucycles();
// 			Hash2(Auth_a, Ke_a, conv<vec_ZZ>(c_b), conv<vec_ZZ>(k_a));
// 				t2 = cpucycles();
// 				ta_h2 = t2-t1;

// 			if (IsZero(Auth_a - Auth_b)){

// 					t1 = cpucycles();
// 				Hash1(sk_a, Ke_a, conv<vec_ZZ>(c_b), Auth_b, conv<vec_ZZ>(k_b));
// 					t2 = cpucycles();
// 					ta_h1 = t2-t1;

// 				//cout <<"SK_a = "<<sk_a <<"\n";

// 			}
// 			else{cout << "Alice: Abort! (H2 Hashes Didn't Match!)"; return; }
			

// 		}
// 		else{ cout << "Alice: Abort!"; return; }	


// 	//CHECK IF EVERYTHING WORKED!
// 		if (IsZero(sk_a - sk_b)){
// 			cout << "\nsk_a == sk_b; Successful AKE!\n";
// 		}


// 		cout << "N=" << N0 <<", q=" << q0 << "\n";
// 		cout << "KEM_Norm=" << KEM_NORM << ", DS_Sigma=" << DS_SIGMA <<"\n\n";
// 		cout << "\n\n========CPU Cycles=======\n";
// 		cout << "\n==ALICE==\n";
// 		cout << "- DS_keygen  :"<< ta_ds_keygen		<< "\n\n";
// 		cout << "- KEM_keygen :"<< ta_kem_keygen	<< "\n";
// 		cout << "- DS_sign    :"<< ta_ds_sig		<< "\n";
// 		cout << "- DS_verify  :"<< ta_ds_ver		<< "\n";
// 		cout << "- KEM_dec    :"<< ta_kem_dec		<< "\n";
// 		cout << "- H2         :"<< ta_h2			<< "\n";
// 		cout << "- H1         :"<< ta_h1			<< "\n";
// 		cout << " TOTAL (ALL) :"<< ta_ds_keygen + ta_kem_keygen + ta_ds_sig + ta_ds_ver + ta_kem_dec + ta_h2 + ta_h1 << "\n";
// 		cout << " TOTAL (AKE) :"<< ta_kem_keygen + ta_ds_sig + ta_ds_ver + ta_kem_dec + ta_h2 + ta_h1 << "\n";

// 		cout << "\n==BOB==\n";
// 		cout << "- DS_keygen  :"<< tb_ds_keygen		<< "\n\n";
// 		cout << "- DS_verify  :"<< tb_ds_ver		<< "\n";
// 		cout << "- KEM_enc    :"<< tb_kem_enc		<< "\n";
// 		cout << "- H2         :"<< tb_h2			<< "\n";
// 		cout << "- DS_sign    :"<< tb_ds_sig		<< "\n";
// 		cout << "- H1         :"<< tb_h1			<< "\n";
// 		cout << " TOTAL (ALL) :"<< tb_ds_keygen + tb_ds_ver + tb_kem_enc + tb_h2 + tb_ds_sig + tb_h1 << "\n";
// 		cout << " TOTAL (AKE) :"<< tb_ds_ver + tb_kem_enc + tb_h2 + tb_ds_sig + tb_h1 << "\n";

// }

int main(){

	

    srand(rdtsc());

    if (false){
	    cout <<"\n\n--RUNNING THE DS/KEM Examples--\n\n";
	    //Test the DS
	    run_DS_example();

	    //Test the KEM
	    run_KEM_example();
	}	

    if (false){
	    //AKE Example
	    cout <<"\n\n--RUNNING THE AKE Examples--\n\n";
	    AKE_timed_exampleMR();
	    //AKE_clocked_example();
	}

	if (false){
		cout <<"\n\n--WORKSPACE: MESSAGE-RECOVERY--\n\n";
		run_DS_example();
	}

	if (true){
		cout << "\n\n--WORKSPACE: GENERATE HUFFMAN--\n\n";
		SampleSave();
	}
    return 0;
}





