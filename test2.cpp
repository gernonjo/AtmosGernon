//-*-mode:c++; mode:font-lock;-*-
/**
 * \class AtmTest
 * \ingroup
 * \brief This code is used to test the modified kasatm/atmosphere classes
 * they have been modifed by GHS to use the CORSIKA AtmProf tables (or the 
 * original US76)
 *
 * Original Author: Glenn H. Sembroski * $Author: cduke $
 * $Date: 2011/08/10 20:03:53 $
 * $Revision: 1.1 $
 * $Tag$
 *
 **/

//Written by:
// G.H.Sembroski
//Physics Dept.
//Purdue Univ.
//West Lafayette, In. 479096
//sembrosk@physics.purdue.edu
//765-494-5172

// 27/12/10



#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <vector>

using namespace std;

#include "TMath.h"
#include "TRandom3.h"
#include "kasatmGrISU.h"
#include "TROOT.h"
#include "stdint.h"

const double NM = 1.0E-09;

int main(int argc, char** argv)

{

  // *********************************************************************
  // Since this is just a test program no need for fancy inputs
  // Outputs can be root TTree::ReadFile()compatable.
  // *********************************************************************

  //V.Crab.US76.50km
  initAtmosphere("atmprof6.dat");
  
  double atmStepsKM=1.;
  double atmTopKM = 10.0;
  int numSteps= (atmTopKM/ atmStepsKM) + 1;
  const  double PlankConstant = 4.136e-15;  //eV
  double lambda1 = 300.;  //300 nm.
  double lambda2 = 600.; //600 nm.
  double altM = 10000.0/2.0;
  double wavelengthM = lambda1* (NM);
  double wavelength2M = lambda2* (NM);
  float altF= altM;
  float wavelengthF = wavelengthM;
  double IndexRef = index_of_ref (&wavelengthF, &altF);
  double n = IndexRef;
  double lambda2inv = (1/(lambda2*(NM)));
  double lambda1inv = (1/(lambda1*(NM)));
  double NumPhotons =0;
  double c = TMath::C();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "lambdaMin" <<" = " <<wavelengthM<< endl;
  cout << "lambdaMax" << " = " <<wavelength2M<< endl;
  cout << "wavelength   altF " << wavelengthF << "   " << altF << endl;
  cout << "IndexofRefraction" <<" = " <<IndexRef<< endl;
  cout <<" TMath::C() "<< TMath::C() << endl;
  cout << endl;
  cout << "NumPhotons " << "    wavelength3F"  << "    eVelocity"
       <<"    Beta1" <<"    theta"<< "        n"<< "     xPosition_ground" <<endl; 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // TRandom3 class declaration
  
  TRandom3 ran;
  ran.SetSeed();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  double alt3M=10000.0/2.0;
  double xPosition_ground;
  double eEnergy = 1000000.0/2.0;
  double theta;
  double eVelocity;
    
  for (int i=0; i<50; i++){
  
    NumPhotons= NumPhotons+1;

    double w = ran.Rndm();

    double wavelength_random = w*(lambda1inv-lambda2inv)+lambda2inv;

    double wavelength3M =(1/wavelength_random);

    float altF=(float)alt3M;

    float wavelength3F = (float)wavelength3M;

    double IndexRef1=index_of_ref (&wavelength3F, &altF);

    double n1 = IndexRef1;

    eVelocity = sqrt(1-(.261121/eEnergy))*c;

    double Beta1 = (eVelocity/c);

    double g2 = (1/(Beta1*n1));

    theta= (acos(g2));

    double sin2 = sin(theta);

    sin2 = (sin2*sin2);

    xPosition_ground=alt3M*tan(theta);
 
       
    cout << NumPhotons << "              " <<  wavelength3F << "    " << eVelocity <<"     "<< Beta1 <<"     "<< theta*TMath::RadToDeg() <<"     " << n << "    "<< xPosition_ground <<  endl;
     }
  return 0; }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

// **************************************************************************


