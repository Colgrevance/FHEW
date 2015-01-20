#ifndef FHEW_H
#define FHEW_H

#include "params.h"
#include "LWE.h"

namespace FHEW {

  typedef Ring_FFT  ct_FFT[K2][2];                        // Ciphertext in FFT form
  typedef ct_FFT* BootstrappingKey[n][BS_base][BS_exp];
  typedef struct {
    BootstrappingKey BSkey;
    LWE::SwitchingKey KSkey;
  } EvalKey;

  void Setup();
  void KeyGen(EvalKey* EK, const LWE::SecretKey LWEsk);
  void HomNAND(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2);
  
  void HomNANDProf(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2, long long * initP, long long * ACCP, long long * KeySwitchP, long long * ModSwitchP, long long * inACC1T, long long * inACC2T, long long * inACC3T, long long * inACC4T, int * numACC);                                          // Profiling variant of HomNAND.

  void AddToACCProf(ct_FFT ACC, ct_FFT C, long long * inACC1, long long * inACC2, long long * inACC3, long long * inACC4);      // Profiling variant of AddToACC

  void fwrite_ek(const EvalKey& EK, FILE* f);
  EvalKey* fread_ek(FILE* f);
}

#endif
