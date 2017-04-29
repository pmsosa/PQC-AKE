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


#include <iostream>
#include <map>

#include "params.h"
#include "Random.h"
#include "Algebra.h"
#include "DigitalSignature.h"


using namespace std;

std::map<int,std::string> enc_huffman;
std::map<std::string,int> dec_huffman;


const int ERROR = q0+100;

//Generate many s_1, and s_2 to then figure out a huffman encoding with python.
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

void LoadHuffman(){
    ifstream file;
    file.open("huffman/huffman.ser");
    if (!file.is_open()) return;

    int a;
    string b;
    while (file >> a){
        file >> b;
        enc_huffman[a] = b;
        dec_huffman[b] = a;
    }
}



//Encode Integer to Huffman code
string Encode(int a){
    try{
        return enc_huffman[a];
    }
    catch(int e){
        cout << "Huffman Error!\n";
        return "";
    }
}


//Decode Huffman code to integer
int Decode(string b){
    try{
        if (b != Encode(0) and dec_huffman[b] == 0){
            return ERROR;
        }
        else{
            return dec_huffman[b];
        }
    }
    catch(int e){
        cout << "Huffman Error!\n";
        return ERROR;   
    }
}

vec_ZZ DecodeFullString(string str){

    vec_ZZ a;
    a.SetLength(1);

    int i = 0;
    unsigned int offset = 0;
    int size = 8;

    string test = "";

    while(offset < str.size()){
        int r = Decode(str.substr(offset,size));
        //cout << str.substr(offset,size) << " = " << r << "\n";
        if (r == ERROR){
            size += 1;
        }
        else{
            offset += size;
            size = 8;
            a[i] = r;
            i +=1;
            a.SetLength(a.length()+1);
            //cout << "i:" << offset << "\n";
            //cout << str.size();

        }
    }
    return a;

}

string EncodeFullVector(vec_ZZ a){

    string result = "";
    for (int i=0; i < a.length(); i++){
        result.append(Encode(conv<int>(a[i])));
        //cout << "      a:"<<a[i] << " =" << Encode(conv<int>(a[i])) <<"\n";
    }
    return result;

}

void testHuffman()
{

    LoadHuffman();
    cout << "Loaded!\n";
    string s1 = Encode(-1);
    cout << "s1: " << s1 << "\n";
    string s2 = Encode(123);
    cout << "s2: " << s2 << "\n";
    string s3 = Encode(0);
    cout << "s3: " << s3 << "\n";
    string s4 = Encode(-412);
    cout << "s4: " << s4 << "\n";

    cout << "Decode('132222'): " << Decode("132222") << "\n";

    string s5 = s1 + s2 + s3 + s4;

    vec_ZZ a;
    a.SetLength(4);
    a[0] = -1;
    a[1] = 123;
    a[2] = 0;
    a[3] = -412;
    cout << "Vect:" << a << "\n";
    cout << "Enc:"<<EncodeFullVector(a) << "\n";
    cout << "Dec:"<<DecodeFullString(EncodeFullVector(a)) << "\n";


    cout << "s5: " << s5 << "\n";
    
    cout << DecodeFullString(s5) << "\n";






  
  // mymap[-1]="1000011";
  // mymap[-2]="00101";
  // mymap[-3]=mymap[-1];

  // std::cout << "mymap[-3] is " << mymap[-3] << '\n';
  // std::cout << "mymap[-2] is " << mymap[-2] << '\n';
  // std::cout << "mymap[-2] is " << mymap[-2] << '\n';
  // std::cout << "mymap[-1] is " << mymap[-1] << '\n';

  // std::cout << "mymap now contains " << mymap.size() << " elements.\n";
}
