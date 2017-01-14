#include<stdio.h>       //Printing 
#include<stdlib.h>      //Printing
#include<math.h>        //Random    
#include<time.h>        //Seed random
#include<complex.h>     //Creal
#include<string.h>      //Memcopy
#include<fftw3.h>       //FTT

//Compile: gcc test.c -lfftw3 -lm

void test1(); //Simple Test were we grab a number get the fft and then we ifft it.
void mult(); //Actual Poly Multiplication.

int main(void)
{
    
    //test1();
    mult();
}

void mult(){
    
    //TODO:
    //1. Cleanup
    //2. Turn into function that takes two ints.
    //3. Make it work with uinit_16 (which is what newhope uses)
    //4. Optimize

    //Setup
    int N = 4;
    double *a, *b, *in;                 // Actual Polynomials (IN)
    fftw_complex *a2, *b2, *out, *c;    // Actual Polynomials (OUT)
    fftw_plan rfft, irfft;              // FFTW Plans
    
    //Space Allocation
    a = fftw_malloc(sizeof(double)*(2*N-1));
    b = fftw_malloc(sizeof(double)*(2*N-1));
    in = fftw_malloc(sizeof(double)*(2*N-1));
    
	a2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	b2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(2*N-1));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(2*N-1));
    
    //Create Plans for rfft and irfft
    // Overhead vs Suboptimanl = FFTW_ESTIMATE | FFTW_MEASURE
	rfft = fftw_plan_dft_r2c_1d(2*N-1, in, out, FFTW_ESTIMATE);
    irfft = fftw_plan_dft_c2r_1d(2*N-1, out, in, FFTW_ESTIMATE);
	

    //Randomly build coefficients
    srand(5); //srand(time(NULL));
    
	for(int i = 0; i < 2*N-1; i++)
	{
		if (i < N){
            a[i] = round(rand()/(RAND_MAX+1.)*10+1);
            b[i] = round(rand()/(RAND_MAX+1.)*10+1);
        }
        else{
            a[i] = 0.;
            b[i] = 0.;
        }
	}
    
    //Print Coefficients
	printf("Polynomial  coefficients:\n");
	for(int i = 0; i < N; i++)
	{
		printf("%2d %11.7f %11.7f\n", i, a[i], b[i]);
	}
    
    //rfft(a)
    memcpy (in,a,sizeof(a)*(2*N-1));
    fftw_execute(rfft);
    memcpy (a2,out,sizeof(fftw_complex)*(2*N-1));

    //rfft(b)
    memcpy (in,b,sizeof(b)*(2*N-1));
    fftw_execute(rfft);
    memcpy (b2,out,sizeof(fftw_complex)*(2*N-1));
    
    // c = irfft(rfft(a')*rfft(b'))
    for(int i = 0; i < 2*N-1; i++)
	{
		out[i] = (a2[i])*(b2[i]);
	}
    
    printf("Temp Check:\n");
    for(int i = 0; i < 2*N-1; i++)
	{
		printf("%2d %11.7f %11.7f %11.7f\n", i, a2[i], b2[i],out[i]);
	}
    
    fftw_execute(irfft);

    //Print Result
	printf("\nPolynomial coefficients of the product:\n");
    for(int i = 0; i < 2*N-1; i++)
	{
		printf("%2d %11.7f\n", i, in[i]/(2*N-1));
	}
}
    
    
    
    

void test1()
{
    //Params & Constant Setup
    int N = 4; //N - Size of Polynomials
    double *in;
	fftw_complex *out;
	fftw_plan rfft, irfft; //FFTW Plan
    
    //Allocating Space for Polynomials
	in = fftw_malloc(sizeof(double)*N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    
    
    //Create Plans for rfft and irfft
    // Overhead vs Suboptimanl = FFTW_ESTIMATE | FFTW_MEASURE
	rfft = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    irfft = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
    
    //Result is unormalized so it is scaled by n.
	

    //Aliases for Accessing Complex Transformation Outputs
    //A = (fftw_complex*) &in[0][0];
    //B = (fftw_complex*) &out[0][0];
    
    //A[0][0] = 0
    //A[1][0] = 1
    //A[2][0] = 0
    //A[3][0] = 0
    
    printf("Loc: %d,%d,%d,%d,\n",&in[0],&in[1],&in[2],&in[3]);
    printf("Data: %f,%f,%f,%f,\n\n",in[0],in[1],in[2],in[3]);
    
    in[0] = 1;
    in[1] = 2;
    in[2] = 3;
    in[3] = 4;
    
    printf("Loc: %d,%d,%d,%d,\n",&in[0],&in[1],&in[2],&in[3]);
    printf("Data: %f,%f,%f,%f,\n\n",in[0],in[1],in[2],in[3]);
    
    
    
    //Multiply    
    fftw_execute(rfft);
    fftw_execute(irfft);
    
    printf("Loc: %d,%d,%d,%d,\n",&in[0],&in[1],&in[2],&in[3]);
    printf("Data: %f,%f,%f,%f,\n\n",in[0],in[1],in[2],in[3]);
	
    
    
    
    
    //Cleanup
    fftw_destroy_plan(rfft);
    fftw_destroy_plan(irfft);
	fftw_free(in);
	fftw_free(out);
}


