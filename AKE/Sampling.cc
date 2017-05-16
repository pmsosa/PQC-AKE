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

#include "Sampling.h"
#include "params.h"

using namespace std;
using namespace NTL;



//==============================================================================
// Takes in input a random value and samples from distribution D_{\sigma_2}^+,   
// Samples an element in Z+ with probability proportionnal to 2^{-x^2}       
//==============================================================================
unsigned int Sample0(unsigned long alea)
{
    if((alea&1UL)==0UL)
    {
        return 0;
    }
    unsigned int i;
    unsigned int k = 1;
    unsigned long mask=0;
    unsigned long aux;
    for(i=1; i<1000;)
    {
        aux = (alea&mask);
        alea = (alea>>k);     
        if(aux)
        {
            return Sample0(alea);
        }
        else
        {
            if((alea&1UL)==0UL)
            {
                return i;
            }
        }
        i++;
        k += 2;
        mask = (mask<<2)|6UL;
    }
    cout << "ERROR" << endl;
    return 999999;
}



//==============================================================================
// Samples from distribution D_{k\sigma_2}^+, ie 
// Samples an element in Z+ with probability proportionnal to 2^{-(x/k)^2} 
//==============================================================================
unsigned int Sample1(const unsigned int k)
{
    unsigned int x, y, z;
    unsigned long alea = rand();

    x = Sample0(alea);
    y = rand()%k;
    z = k*x + y;
    RR_t w = y*( (z<<1) - y );
    RR_t borne =  LDRMX / exp( w*log_2/(k*k) );
    alea = rand();
    if(alea>borne)
    {
        return Sample1(k);
    }
    else
    {
        return z;
    }
    cout << "ERROR" << endl;
    return 999999;
}




//==============================================================================
// Samples from distribution D_{k\sigma_2}, ie                
// Samples an element in Z with probability proportionnal to 2^{-(x/k)^2} 
//==============================================================================
signed int Sample2(const unsigned int k)
{
    signed int signe;
    signed int x;
    unsigned long alea = rand();
    while(1)
    {
        x = Sample1(k);
        if( (x!=0) || ((alea&1)==1) )
        {
            alea >>= 1;
            signe = 1 - 2*(alea&1);
            x *= signe;
            return x;
        }
        alea >>= 1;
    }
}


//==============================================================================
// Samples from distribution D_{sigma}, ie                                       
// Samples an element in Z with probability proportionnal to e^{-x^2/2*(sigma^2)}
//==============================================================================
signed int Sample3(const RR_t sigma128)
{
    signed int x;
    double alea, borne;

    const RR_t sigma = sigma128;
    const unsigned long k = ( (unsigned long) ceil( (RR_t) sigma/sigma_1 ) );
    while(1)

    {
        x = Sample2(k);
        alea = ((RR_t)rand()) / LDRMX;
        borne = exp( -x*x*( 1/(2*sigma*sigma) - 1/(2*k*k*sigma_1*sigma_1) )   );
        assert(borne<=1);
        if(alea<borne)
        {
            return x;
        }
    }
}


//==============================================================================
// Samples from distribution D_{c,sigma}, ie                                              
// Samples an element in Z with probability proportionnal to e^{-(c-x)^2/2*(sigma^2)}    
//==============================================================================
signed int Sample4(RR_t c, RR_t sigma)
{
    RR_t alea, borne;
    signed int x;
    unsigned int coin;

    const signed int intc = ( (signed int) floor(c) );
    const RR_t fracc = c-intc;
    coin = rand();
    const RR_t denom = 1/(2*sigma*sigma);

    while(1)
    {
        x = Sample3(sigma);
        x += (coin&1);
        if(abs(x)>8){cout << x << endl;}
        coin >>= 1;
        borne = exp(-(x-fracc)*(x-fracc)*denom)/ ( exp(-x*x*denom) + exp(-(x-1)*(x-1)*denom) );

        assert(borne<1);
        alea = ( (RR_t)rand() ) / LDRMX;
        if(alea<borne)
        {
            return (x+intc);
        }
    }
}
