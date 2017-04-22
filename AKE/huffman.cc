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
#include <iostream>
#include <fstream>
#include <string>

#include "Sampling.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"
#include "DigitalSignature.h"


using namespace std;

void SampleSave(){

    ZZX Ks[2], s[2];;
    ZZ_pX Kv;
    vec_ZZ r, m2;
    MSK_Data * MSKD = new MSK_Data;

    SigKeyGen(Ks,Kv,MSKD);

    /// RUNNING THE SIGN ALGORITHM
    int i = 0;
    int write_limit = 140967296;
    int written = 0;
    bool open_new = true;
    int filelimit = 14;
    int progress = 0;
    ofstream myfile;

    //
    while(i <= filelimit){
        //Open File
        if (open_new){
            myfile.open("/media/sf_CS276/example"+std::to_string(i)+".txt");
            open_new = false;
        }

        //Generate s_1 and s_2
        vec_ZZ msg = RandomVector();
        Sign2(s,m2,msg,MSKD);

        //Save s_1 and s_2 in file
        for (int j=0; j < deg(s[0]);j++){
            myfile << conv<int>(s[0][j]) << "\n";
            written += 1;
            progress += 1;
        }
        for (int j=0; j < deg(s[1]);j++){
            myfile << conv<int>(s[1][j]) << "\n";
            written += 1;
            progress += 1;
        }

        if (progress >= write_limit*0.001){
            progress = 0;
            int x = written*100/write_limit;
            cout  << x << " ";
        }

        // if ((written >= write_limit*0.1*progress)==0){
        //     cout << written/write_limit << " ";
        //     progress += 1;
        // }

        if (written >= write_limit){
            myfile.close();
            i += 1;
            written = 0;
            open_new = true;
            cout << "\nFinished Writing: " << i-1 << "\n\n\n";
        }

    }

}