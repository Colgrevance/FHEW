#include <iostream>
#include <cstdlib>
#include "LWE.h"
#include "FHEW.h"
#include "distrib.h"

#include "time_profiler.h"

using namespace std;

void help(char* cmd) {
  cerr << "\nusage: " << cmd << " n\n\n" 
  << "  Run it again. Do not forget to write a number after the command. \n\n";
  exit(0);
}

int main(int argc, char *argv[]) {
  if (argc != 2) help(argv[0]);
  int count = atoi(argv[1]); 

  cerr << "Setting up FHEW \n";
  FHEW::Setup();

  cerr << "Generating secret key ... ";
  LWE::SecretKey LWEsk;
  LWE::KeyGen(LWEsk);
  cerr << " Done.\n";
  
  cerr << "Generating evaluation key ... this may take a while ... ";
  FHEW::EvalKey EK;
  FHEW::KeyGen(&EK, LWEsk);
  cerr << " Done.\n\n";
  
  cerr << "Testing homomorphic NAND " << count << " times.\n"; 
  cerr << "Circuit shape : (a NAND b) NAND (c NAND d)\n\n";

  int v1,v2;  
  LWE::CipherText e1, e2, e12;


  
  long long initAv = 0, ACCAv = 0, KeySwitchAv = 0, ModSwitchAv = 0, inACC1T = 0, inACC2T = 0, inACC3T = 0, inACC4T = 0;     // Profile values for the average of the functions.

  for (int i = 0; i < count; ++i) {

    v1 = rand()%2;  
    v2 = rand()%2;
    LWE::Encrypt(&e1, LWEsk, v1);
    LWE::Encrypt(&e2, LWEsk, v2);
    
    cerr << "Enc(" << v1 << ")  NAND  Enc(" << v2 << ")  =  ";


    long long initP, ACCP, KeySwitchP, ModSwitchP, inACC1P, inACC2P, inACC3P, inACC4P; // Profile values for the functions inside NAND.
    int numACC;
    FHEW::HomNANDProf(&e12, EK, e1, e2, &initP, &ACCP, &KeySwitchP, &ModSwitchP, &inACC1P, &inACC2P, &inACC3P, &inACC4P, &numACC);                                     // call the profiling variant of HomNAND
    int v12 = LWE::Decrypt(LWEsk, e12);

    cerr << "Enc(" << v12 << ")";
    cerr << endl;

    if (1 - v1*v2 != v12) { 
      cerr << "ERROR at iteration " << i+1 << "\n"; 
    exit(1);
    }

    cout << "The profile for the " << i+1 << "th iteration of HomNAND is: \n";
    cout << "For InitializeACC:      " << initP << "\n";
    cout << "For AddToACC (overall): " << ACCP << ", divided into: \n";
    cout << "                        " << inACC1P << " for the first loop\n";
    cout << "                        " << inACC2P << " for the second loop\n";
    cout << "                        " << inACC3P << " for the third loop\n";
    cout << "                        " << inACC4P << " for the fourth loop\n";
    cout << "For KeySwitch:          " << KeySwitchP << "\n";
    cout << "For ModSwitch:          " << ModSwitchP << "\n";
    cout << "The number of ACC is " << numACC << "\n"; 
    cout << endl;


    initAv = initAv + initP;
    ACCAv = ACCAv + ACCP;
    KeySwitchAv = KeySwitchAv + KeySwitchP;
    ModSwitchAv = ModSwitchAv + ModSwitchP;
    inACC1T = inACC1T + inACC1P;
    inACC2T = inACC2T + inACC2P;
    inACC3T = inACC3T + inACC3P;
    inACC4T = inACC4T + inACC4P;



  }

  cerr << "\nPassed all tests!\n\n";

  cout << "The average for InitializeACC is " << initAv / count << "\n";
  cout << "The average for AddToACC is " <<  ACCAv / count << "\n";
  cout << "The average for KeySwitch is " << KeySwitchAv / count << "\n";
  cout << "The average for ModSwitch is " << ModSwitchAv / count << "\n\n";
  cout << "The average for the first ACC loop is " << inACC1T / count << "\n";
  cout << "The average for the second ACC loop is " << inACC2T / count << "\n";
  cout << "The average for the third ACC loop is " << inACC3T / count << "\n";
  cout << "The average for the fourth ACC loop is " << inACC4T / count << "\n\n\n";


}