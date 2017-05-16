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

#include "FFT.h"
#include "params.h"

using namespace std;
using namespace NTL;


//FFT Multiply
void FFTmultiply(ZZX& result, ZZX& a, ZZX& b){
    unsigned long i;
    CC_t a_FFT[N0], b_FFT[N0], result_FFT[N0];
    ZZXToFFT(a_FFT, a);
    ZZXToFFT(b_FFT, b);

    for(i=0; i<N0; i++)
    {
        result_FFT[i] = a_FFT[i]*b_FFT[i];
    }

    FFTToZZX(result,result_FFT);

}


void FFTStep(CC_t * const f_fft, RR_t const * const f, const unsigned long N, const CC_t w0)
{
    if(N==1)
    {
        f_fft[0] = f[0];
    }
    else
    {
        if(N==2)
        {
            f_fft[0] = f[0] + ii*f[1];
            f_fft[1] = f[0] - ii*f[1];
        }
        else
        {
            assert(N%2==0);
            RR_t f0[N0/2], f1[N0/2];
            CC_t f0_fft[N0/2], f1_fft[N0/2], wk, w02;
            unsigned int k;
            for(k=0; k<(N/2); k++)
            {
                f0[k] = f[2*k];
                f1[k] = f[2*k+1];
            }

            w02 = w0*w0;
            wk = w0;
            FFTStep(f0_fft, f0, N/2, w02);
            FFTStep(f1_fft, f1, N/2, w02);
            for(k=0; k<N; k++)
            {
                f_fft[k] = f0_fft[k%(N/2)] + wk*f1_fft[k%(N/2)];
                wk *= w02;
            }
        }
    }
}


void ReverseFFTStep(CC_t * const f, CC_t const * const f_fft, const unsigned long N, const CC_t w0)
{
    if(N!=2)
    {
        assert(N%2==0);
        CC_t f0[N0/2], f1[N0/2];
        CC_t f0_fft[N0/2], f1_fft[N0/2], w02, wk;
        unsigned int k;

        w02 = w0*w0;
        wk = w0;

        for(k=0; k<N/2; k++)
        {
            f0_fft[k] = (f_fft[k] + f_fft[k+(N/2)])*0.5l;
            f1_fft[k] = wk*(f_fft[k] - f_fft[k+(N/2)])*0.5l;
            wk *= w02;
        }
        ReverseFFTStep(f0, f0_fft, (N/2), w02);
        ReverseFFTStep(f1, f1_fft, (N/2), w02);

        for(k=0; k<N/2; k++)
        {
            f[2*k] = f0[k];
            f[2*k+1] = f1[k];
        }
    }
    else
    {
        f[0] = (f_fft[0] + f_fft[1])*0.5l;
        f[1] = (f_fft[0] - f_fft[1])*(-0.5l*ii);
    }
}


void MyRealReverseFFT(double * const f, CC_t const * const f_fft)
{
    CC_t fprime[N0];
    unsigned int i;
    ReverseFFTStep(fprime, f_fft, N0, omega_1);

    for(i=0; i<N0; i++)
    {
        f[i] = (fprime[i]).real();
    }
}


void MyIntReverseFFT(long int * const f, CC_t const * const f_fft)
{
    CC_t fprime[N0];
    unsigned int i;

    ReverseFFTStep(fprime, f_fft, N0, omega_1);
    for(i=0; i<N0; i++)
    {
        f[i] = ((long int) round( fprime[i].real() ) );
    }
}


void MyIntFFT(CC_t * f_FFT, const long int * const f)
{
    RR_t f_double[N0];
    unsigned int i;

    for(i=0; i<N0; i++)
    {
        f_double[i] = ( RR_t ( f[i] ) );
    }

    FFTStep(f_FFT, f_double, N0, omega);
}


void ZZXToFFT(CC_t * f_FFT, const ZZX f)
{
    RR_t f_double[N0];
    unsigned int i;

    //assert(deg(f)==N0-1);
    assert(MaxBits(f)<900);
    for(i=0; i<N0; i++)
    {
        if (i <= deg(f)){
            f_double[i] = ( RR_t ( conv<double>(f[i]) ) );
        }
        else{
            f_double[i] = ( RR_t ( conv<double>(0) ) );
        }
    }

    FFTStep(f_FFT, f_double, N0, omega);
}


void FFTToZZX(ZZX& f, CC_t const * const f_FFT)
{
    double f_double[N0];
    unsigned int i;
    MyRealReverseFFT(f_double, f_FFT);

    f.SetLength(N0);
    for(i=0; i<N0; i++)
    {
        f[i] = conv<ZZ>(round(f_double[i]));
    }

}
