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



#include "stdint.h"
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;
#include "TROOT.h"
#include "TVector2.h"
#include "TClass.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TRandom1.h"
#include "kasatmGrISU.h"

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
  double wavelength = 450.;  //450 nm.
  double altM = 100.0;
  double wavelengthM = wavelength* (NM);
  cout << "wavelength" <<" = " <<wavelengthM<< endl;
  float altF= altM;
  float wavelengthF = wavelengthM;
  double IndexRef = index_of_ref (&wavelengthF, &altF);
  cout << "IndexofRefraction" <<" = " <<IndexRef<< endl;
  double pi = 3.14159;
  cout <<" TMath::C() "<< TMath::C() << endl;
  cout << endl;
  cout << " NumPhotons " << "  lambda"  << "     eEnergy"
       <<"     velocitynew" <<"   theta"<< endl; 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // constants
 
  double h = PlankConstant*1e-6;

  double n = IndexRef;
  double lambdaMax = (7e-7); // 700nm
  double lambdaMin = (4e-7); //400nm
  double a = (1.0/137.0);
  double seglength = 50.0;
  double x1;

  ///////////////////////////////////////////////////////////////////

  // velocity


  double c = TMath::C();
  double lambda = lambdaMin; //400nm
  double eEnergy;
  double theta;
  double velocitynew;
  double NumPhotons = 0;


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  


   
  // For loop

  for (int i=0; i<10; i++){
    
    
    eEnergy = eEnergy + 100.0;
    int NumPhotons;
    double eVelocity;
    double xPosition_ground;
    eVelocity = sqrt(1-(.261121/eEnergy))*c;
    double Beta = (eVelocity/c);
    double g = (1/(Beta*n));

    double w;
    double g2 = 0.0;
    theta = 0.0;
    NumPhotons = 0;
    xPosition_ground = 0.0;
    if (eVelocity >= (c/n)){

      g2= g*g;
      double j= (2*pi*a*seglength);
      w= (1/lambdaMin)-(1/lambdaMax);
      theta = (acos(g));
      NumPhotons = j*w*(1-g2);
      if (NumPhotons >0){
        xPosition_ground = altM*tan(theta);
      }
    }
   
    cout << eEnergy << "     " << eVelocity << "    " 
         << theta*TMath::RadToDeg() << "          " << NumPhotons <<"            "
     <<xPosition_ground << endl;
  }
    
  cout << endl; 

  return 0;
}

  



  

// **************************************************************************


