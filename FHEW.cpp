#include <iostream>
#include <cstdlib>
#include "FHEW.h"
#include "FFT.h"
#include <cassert>

#include "time_profiler.h"

using namespace std;

namespace FHEW {

  typedef Ring_ModQ ct_ModQ[K2][2];                       // Ciphertext in coefficient form
  typedef Ring_ModQ dct_ModQ[K2][K2];                     // Decomposed Ciphertext in coeff form
  typedef Ring_FFT  dct_FFT[K2][K2];                      // Decomposed Ciphertext in FFT form
  
  Ring_FFT t_TestMSB;
  

  
  void Setup() {
    
    FFTsetup();                                           // Call FFTsetup. The forward and backward plans are declared at this step.
    Ring_ModQ tmsb;
    tmsb[0]=-1;
    
    for (int i = 1; i < N; ++i) {
      tmsb[i]=1;
    }
    
    FFTforward(t_TestMSB, tmsb);
  }


  void FHEWencrypt(ct_FFT ct, Ring_FFT sk_FFT, int m) { // Used by BootstrappingKeyGen

    Ring_FFT ai;
    ct_ModQ res;
    int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    
    if (mm >= N) { 
      mm -= N; sign = -1; 
    }

    for (int i = 0; i < K2; ++i) {
      
      for (int k = 0; k < N; ++k) {
        res[i][0][k] = rand(); // % Q
      }
      FFTforward(ai, res[i][0]);
      
      for (int k = 0; k < N2; ++k) {
        ai[k] = ((double complex) ai[k]) * ((double complex) sk_FFT[k]);
      }
      FFTbackward(res[i][1], ai);
      
      for (int k = 0; k < N; ++k){
        res[i][1][k] += Sample(Chi1);                   // Add error [a,as+e]
      }

    }


    for (int i = 0; i < K; ++i) {
      res[2*i  ][0][mm] += sign*vgprime[i];             // Add G Multiple
      res[2*i+1][1][mm] += sign*vgprime[i];             // [a,as+e] + X^m *G
    }

    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        FFTforward(ct[i][j], res[i][j]);
      }
    }

  }  

  /*
    \\\\\\\\\\\\\\\\\---- Eval-Key Generation ----///////////////////
  */

  void KeyGen(EvalKey* EK, const LWE::SecretKey LWEsk) {

    LWE::SecretKeyN FHEWsk;
    LWE::KeyGenN(FHEWsk);
    LWE::SwitchingKeyGen(EK->KSkey, LWEsk, FHEWsk);

    Ring_FFT FHEWskFFT;
    FFTforward(FHEWskFFT, FHEWsk);

    for (int i = 0; i < n; ++i) {
      for (int j = 1; j < BS_base; ++j) {
        for (int k = 0; k < BS_exp; ++k) {
          EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
          FHEWencrypt( (*EK->BSkey[i][j][k]), FHEWskFFT, LWEsk[i] * j * BS_table[k] );
        }
      }
    }
  
  }



  void fwrite_ek(const EvalKey& EK, FILE* f) {
    /* --- Write bootstrapping key --- */
    for (int i = 0; i < n; ++i) {
      for (int j = 1; j < BS_base; ++j) {
        for (int k = 0; k < BS_exp; ++k) {
          assert(fwrite(EK.BSkey[i][j][k], sizeof(ct_FFT), 1, f));
        }
      }
    }
    
    /* --- Write switching key --- */
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < KS_base; ++j) {
        for (int k = 0; k < KS_exp; ++k) {
          assert(fwrite(EK.KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
        }
      }
    }

  }




  EvalKey* fread_ek(FILE* f) {

    EvalKey* EK = new EvalKey;
    
    // Read bootstrapping key
    for (int i = 0; i < n; ++i) {
      for (int j = 1; j < BS_base; ++j) {
        for (int k = 0; k < BS_exp; ++k) {
          EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
          assert(fread(EK->BSkey[i][j][k], sizeof(ct_FFT), 1, f));
        }
      }
    }
    cerr << "BSKey Read. \n";
    
    // Read switching key
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < KS_base; ++j) {
        for (int k = 0; k < KS_exp; ++k) {
          EK->KSkey[i][j][k] = new LWE::CipherTextQ;
          assert(fread(EK->KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
        }
      }
    }

    return EK;
  }
  
  
  



  void AddToACC(ct_FFT ACC, ct_FFT C) {
    
    ct_ModQ ct;
    dct_ModQ dct;
    dct_FFT dctFFT;

    // Decompose_ct(dctFFT, ACC);
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        FFTbackward(ct[i][j], ACC[i][j]);
      }
    }


    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j){
        for (int k = 0; k < N; ++k) {
          ZmodQ t = ct[i][j][k] * v_inverse;
          for (int l = 0; l < K; ++l) {
            ZmodQ r = (t << g_bits_32[l]) >> g_bits_32[l];
            //if ((l==2) && (k ==0 ))
              //cout << r << ",";                       // Uncomment these two lines for debugging purposes
            t = (t-r) >> g_bits[l];
            dct[i][j+2*l][k] = r;
          }
        }
      }
    }
    
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < K2; ++j){
        FFTforward(dctFFT[i][j], dct[i][j]);
      }
    }

    // Mult_dct_ct(ACC, dct, C);
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < N2; ++k) {
          ACC[i][j][k] = (double complex) 0.0;
          for (int l = 0; l < K2; ++l) {
            ACC[i][j][k] += ((double complex) dctFFT[i][l][k]) * ((double complex) C[l][j][k]);
          }
        }
      }
    }

  }

  /*
    \\\\\\\\\\\\\\\\\\\\\\\\ ---- AddToACCProf ---- ///////////////////////
    
    The following function is a variant of AddToACC modified just to measure the running time of the four loops inside AddToACC.
  */

  void AddToACCProf(ct_FFT ACC, ct_FFT C, long long * inACC1, long long * inACC2, long long * inACC3, long long * inACC4) {
    
    ct_ModQ ct;
    dct_ModQ dct;
    dct_FFT dctFFT;

    long long t1, t2;                                   // rdtsc-time variables
    //omp_set_num_threads(4);


    t1 = profiler::rdtsc();
    //#pragma omp parallel for num_threads(8) // Apparently this is going to make the thing faster...
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        FFTbackward(ct[i][j], ACC[i][j]);
      }
    }
    t2 = profiler::rdtsc();
    *inACC1 = t2 - t1;


    t1 = profiler::rdtsc();
    //#pragma omp parallel for num_threads(8)
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < N; ++k) {
          ZmodQ t = ct[i][j][k] * v_inverse;
          for (int l = 0; l < K; ++l) {
            ZmodQ r = (t << g_bits_32[l]) >> g_bits_32[l];
            //if ((l==2) && (k ==0 ))
              //cout << r << ",";                       // Uncomment these two lines for debugging purposes
            t = (t-r) >> g_bits[l];
            dct[i][j+2*l][k] = r;
          }
        }
      }
    }
    t2 = profiler::rdtsc();
    *inACC2 = t2 - t1;

    t1 = profiler::rdtsc();
    //#pragma omp parallel for num_threads(8)
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < K2; ++j){
        FFTforward(dctFFT[i][j], dct[i][j]);
      }
    }
    t2 = profiler::rdtsc();
    *inACC3 = t2 - t1;

    t1 = profiler::rdtsc();
    // Mult_dct_ct(ACC, dct, C);
    //#pragma omp parallel for num_threads(8)
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < N2; ++k) {
          ACC[i][j][k] = (double complex) 0.0;
          for (int l = 0; l < K2; ++l) {
            ACC[i][j][k] += ((double complex) dctFFT[i][l][k]) * ((double complex) C[l][j][k]);
          }
        }
      }
    }
    t2 = profiler::rdtsc();
    *inACC4 = t2 - t1;

  }




  
  void InitializeACC(ct_FFT ACC, int m) {	                  // Set a ciphertext to X^m * G (encryption of m without errors)
    
    ct_ModQ res;
    int mm = (((m % q) + q) % q) * (2*N/q);                 // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    if (mm >= N) { mm -= N; sign = -1; }
    
    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < N; ++k) {
          res[i][j][k]=0;
        }
      }
    }

    for (int i = 0; i < K; ++i) {
      res[2*i  ][0][mm] += sign*vgprime[i];                 // Add G Multiple
      res[2*i+1][1][mm] += sign*vgprime[i];                 // [a,as+e] + X^m *G
    }

    for (int i = 0; i < K2; ++i) {
      for (int j = 0; j < 2; ++j) {
        FFTforward(ACC[i][j], res[i][j]);
      }
    }

  }


  
  
  LWE::CipherTextQN* MemberTest(Ring_FFT t, ct_FFT C) {
    
    Ring_FFT temp;
    Ring_ModQ temp_ModQ;
    LWE::CipherTextQN* ct = new LWE::CipherTextQN;
    
    for (int i = 0; i < N2; ++i) { // Compute t*a
      temp[i] = conj(((double complex) C[1][0][i]) * ((double complex)t[i]));  
    }
    
    FFTbackward(temp_ModQ, temp);
    
    for (int i = 0; i < N; ++i) {
      ct->a[i] = temp_ModQ[i];
    }

    for (int i = 0; i < N2; ++i){
      temp[i] = ((double complex) C[1][1][i]) * ((double complex) t[i]);
    }

    FFTbackward(temp_ModQ, temp);
    ct->b = v+temp_ModQ[0];	
    
    return ct;    
  }
  

  /*
    \\\\\\\\\\\\\\\\\\\\\\\\\\ ---- HomNAND ---- /////////////////////////
    
    The following function computes the NAND gate of two bit homomorphically, that is, using the (LWE) ciphertexts corresponding to two plaintexts, it computes the ciphertext corresponding to the NAND of these plaintexts.
  */


  void HomNAND(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2) {

    LWE::CipherText e12;
    
    for (int i = 0; i < n; ++i){
      e12.a[i] = (2*q - (ct1.a[i] + ct2.a[i])) % q;
    }
    e12.b  =  (13 * q / 8) - (ct1.b + ct2.b) % q;

    ct_FFT ACC;
    
    InitializeACC(ACC, (e12.b + q/4) % q);
    

    for (int i = 0; i < n; ++i) {
      int a = (q - e12.a[i] % q) % q;
      
      for (int k = 0; k < BS_exp; ++k, 	a /= BS_base) {
        int a0 = a % BS_base;
        
        if (a0) {
          AddToACC(ACC, *(EK.BSkey[i][a0][k]) );	
        
        }
      }
    }

    LWE::CipherTextQN *eQN = MemberTest(t_TestMSB, ACC);
    LWE::CipherTextQ eQ;
    LWE::KeySwitch(&eQ, EK.KSkey, *eQN);
    LWE::ModSwitch(res, eQ);
  }


  /*
    \\\\\\\\\\\\\\\\\\\\\\\\ ---- HomNANDProf ---- ///////////////////////
    
    The following function is a variant of HomNAND modified just to measure the running time of the different functions that HomNAND is using.
  */

  void HomNANDProf(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2, long long* initP, long long* ACCP, long long* KeySwitchP, long long* ModSwitchP, long long * inACC1P, long long * inACC2P, long long * inACC3P, long long * inACC4P, int * numACC) {

    /* Variables created for profiling purposes */
    unsigned long long t1, t2;
    long long inACC1temp = 0, inACC2temp = 0, inACC3temp = 0, inACC4temp = 0, inACC1, inACC2, inACC3, inACC4;
    int nac = 0;


    LWE::CipherText e12;
    
    for (int i = 0; i < n; ++i){
      e12.a[i] = (2*q - (ct1.a[i] + ct2.a[i])) % q;
    }
    e12.b  =  (13 * q / 8) - (ct1.b + ct2.b) % q;

    ct_FFT ACC;

    /*
      Profiling InitializeACC
    */

    t1 = profiler::rdtsc();
    InitializeACC(ACC, (e12.b + q/4) % q);
    t2 = profiler::rdtsc();
    *initP = t2 - t1;


    /*
      Profiling ACC
    */
    
    t1 = profiler::rdtsc();
    for (int i = 0; i < n; ++i) {
      int a = (q - e12.a[i] % q) % q;
      
      for (int k = 0; k < BS_exp; ++k,  a /= BS_base) {
        int a0 = a % BS_base;
        
        if (a0) {
          AddToACCProf(ACC, *(EK.BSkey[i][a0][k]), &inACC1, &inACC2, &inACC3, &inACC4);  
            inACC1temp = inACC1temp + inACC1;
            inACC2temp = inACC2temp + inACC2;
            inACC3temp = inACC3temp + inACC3;
            inACC4temp = inACC4temp + inACC4;
            nac = nac + 1;
        }
      }
    }
    t2 = profiler::rdtsc();
    *ACCP = t2 - t1;

    *inACC1P = inACC1temp;        // Pass the running time of the first loop
    *inACC2P = inACC2temp;        // Pass the running time of the second loop
    *inACC3P = inACC3temp;        // Pass the running time of the third loop
    *inACC4P = inACC4temp;        // Pass the running time of the fourth loop
    *numACC = nac;                // Pass the number of AddToACC calls


    /*
      Profiling KeySwitch
    */

    LWE::CipherTextQN *eQN = MemberTest(t_TestMSB, ACC);
    LWE::CipherTextQ eQ;

    t1 = profiler::rdtsc();
    LWE::KeySwitch(&eQ, EK.KSkey, *eQN);
    t2 = profiler::rdtsc();
    *KeySwitchP = t2 - t1;


    /*
      Profiling ModSwitch
    */

    t1 = profiler::rdtsc();
    LWE::ModSwitch(res, eQ);
    t2 = profiler::rdtsc();
    *ModSwitchP = t2 - t1;

  }

  
}
