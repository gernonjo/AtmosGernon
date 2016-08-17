

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
  //**************************************************************************
  ///////////////////////////////////////////////////////////////////////////
  // This is my program for the position of cherenkov photons on the ground.
  //**************************************************************************  

  initAtmosphere("atmprof6.dat");
 
  double atmStepsKM=1.;
  double atmTopKM = 10.0;
  int numSteps= (atmTopKM/ atmStepsKM) + 1;
  const  double PlankConstant = 4.136e-15;  //eV
  double lambda1 = 300.;  //300 nm.
  double lambda2 = 600.; //600 nm.
  double altM = 10000.0;
  double wavelengthM = lambda1* (NM);
  double wavelength2M = lambda2* (NM);
  float altF= altM;
  float wavelengthF = wavelengthM;
  double IndexRef = index_of_ref (&wavelengthF, &altF);
  double n = IndexRef;
  double lambda2inv = (1/(lambda2*(NM)));
  double lambda1inv = (1/(lambda1*(NM)));
  double c = TMath::C();
  double pi = TMath::Pi();
  double seglength= 50;
  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout << "lambdaMin" <<" = " <<wavelengthM<< endl;
  cout << "lambdaMax" << " = " <<wavelength2M<< endl;
  cout << "wavelength   altF " << wavelengthF << "   " << altF << endl;
  cout << "IndexofRefraction" <<" = " <<IndexRef<< endl;
  cout <<"TMath::C() "<< TMath::C() << endl;
  cout <<"TMath::Pi() "<< TMath::Pi() << endl;
  cout << endl;
  cout << "wavelegth3F " << "    photalt"  << "    eVelocity"
       <<"      NumPhotons" <<"     theta"<< "          n"<< "      xPosition_ground" <<endl; 

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 // TRandom3 class declaration
  
  TRandom3 ran;
  ran.SetSeed();
  
  // TVector3 declaration
  TVector3 v;
  TVector3 v1;
  TVector3 v2;
  TVector3 v3;
 
   

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double NumPhotons;
  double SetTheta;
  double alt3M = 10000.0;
  double xPosition_ground;
  double eEnergy = 1000000.0/2.0;
  double theta;
  double eVelocity;
  double a = (1.0/137.0);
  double j;
  double lamb;
  double mag = v.Mag()+1;
 
  double thet = v.Theta();
  
 
    
    double w = ran.Rndm();
    double wavelength_random = w*(lambda1inv-lambda2inv)+lambda2inv;
    double wavelength3M =(1/wavelength_random);
    float altF=(float)alt3M;
    float wavelength3F = (float)wavelength3M;
    double IndexRef1=index_of_ref (&wavelength3F, &altF);
    double n1 = IndexRef1;
  for (int i=0; i<50; i++){

    eVelocity = sqrt(1-(.261121/eEnergy))*c;
    double Beta1 = (eVelocity/c);
    double g = (1/(Beta1*n1));
    double g2 = 0.0;
    double photalt = w*(seglength)+(alt3M-seglength);
    double azim = w*360;
    double z;
    double x;
    double y;
   

    if (eVelocity >= (c/n1)){

 
      g2=g*g;
      j= (2*pi*a*seglength);
      lamb = (1/wavelengthM)-(1/wavelength2M);
      NumPhotons = j*lamb*(1-g2);
      theta= (acos(g2));
      double sin2 = sin(theta);
      sin2 = (sin2*sin2);
      
      x = sin(theta)*cos(azim);
      y = sin(theta)*sin(azim);
      z = cos(theta);
      //  v1.XYZ(x,y,z);

      //   cout<< x<< "  " << y << "  " << z << "  " << endl;
     
  

      xPosition_ground = 0.0;
      if (NumPhotons >0){

        xPosition_ground=alt3M*tan(theta);
    }
    }
    
    cout << wavelength3F << "      " <<  photalt << "    " << eVelocity <<"     "<< NumPhotons <<"     "<< theta*TMath::RadToDeg() <<"     " << n << "    "<< xPosition_ground <<  endl; 
  }
  return 0; }
