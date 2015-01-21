/*
  ////////////////////---- FFT Constructor ---- \\\\\\\\\\\\\\\\\\\\

    
*/

 #include <iostream>                                     // Uncomment this line for debugging purposes
#include "time_profiler.h"

#include <complex.h>
#include <fftw3.h>
#include "FFT.h"
#include "params.h"



double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;
  
void FFTsetup() {
  in = (double*) fftw_malloc(sizeof(double) * 2*N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N + 2));
  /*
    The following two lines form the planner
    The last entry is the planning-rigor flag. The possible values are 
      -FFTW_ESTIMATE
      -FFTW_MEASURE
      -FFTW_PATIENT
      -FFTW_EXHAUSTIVE
  */
  plan_fft_forw = fftw_plan_dft_r2c_1d(2*N, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2*N, out, in,  FFTW_PATIENT);
}

/*  
  The variable in ends up of the form:
    [val[0], val[1],...,val[N], 0.0, 0.0,...,0.0]
*/
  
void FFTforward(Ring_FFT res, const Ring_ModQ val) {
  for (int k = 0; k < N; ++k)	{
    in[k] = (double) (val[k]);
    in[k+N] = 0.0;			
  }
  
  fftw_execute(plan_fft_forw); 
  for (int k = 0; k < N2; ++k) {                        // Transforms the output of the FFT into an array of long complexes
    res[k] = (double complex) out[2*k+1];
  }

  //for (int i = 0; i<N2; ++i) std::cout << res[i];     // Uncomment this line for debugging purposes

  /*
    Uncomment the following lines to free the pointers
    TODO: Test this.
  */

  //fftw_destroy_plan(plan_fft_forw); 
  //fftw_free(in); fftw_free(out);
}

void FFTbackward(Ring_ModQ res, const Ring_FFT val){
  
  //for ( int i = 0; i<N; ++i) std::cout << res[i];      std::cout << "\n";    // Uncomment this line for debugging purposes

  for (int k = 0; k < N2; ++k) {
    out[2*k+1] = (double complex) val[k]/N;
    out[2*k]   = (double complex) 0;
  }

  fftw_execute(plan_fft_back); 
  for (int k = 0; k < N; ++k)	{                         // Transforms the output of the reverse FFT into an aray of long integers
    res[k] = (long int) round(in[k]);
  }

  //for ( int i = 0; i<N; ++i) std::cout << res[i];       std::cout << "\n";     // Uncomment this line for debugging purposes

  /*
    Uncomment the following lines to free the pointers
    TODO: Test this.
  */

  //fftw_destroy_plan(plan_fft_forw); 
  //fftw_free(in); fftw_free(out);
}
