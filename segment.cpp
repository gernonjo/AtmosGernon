

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
#include "TRotation.h"

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
  double seglength= 10000;
  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << endl;
  cout << "lambdaMin" <<" = " <<wavelengthM<< endl;
  cout << "lambdaMax" << " = " <<wavelength2M<< endl;
  cout << "wavelength   altF " << wavelengthF << "   " << altF << endl;
  cout << "IndexofRefraction" <<" = " <<IndexRef<< endl;
  cout <<"TMath::C() "<< TMath::C() << endl;
  cout <<"TMath::Pi() "<< TMath::Pi() << endl;
  cout << endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 // TRandom3 class declaration
  TRandom3 ran;
  ran.SetSeed();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double NumPhotons;
  double alt3M = 10000.0;
  double eEnergy = 1000000.0/2.0;
  double eVelocity;
  double a = (1.0/137.0);
  double j;
  double lamb;
  double theta;
  

  eVelocity = sqrt(1-(.261121/eEnergy))*c;
  double Beta1 = (eVelocity/c);
  double g = (1/(Beta1*n));
  double g2 = g*g;
  theta = acos(g2); 
  cout<< theta << endl;

  if (eVelocity >= (c/n)){

    g2=g*g;
    j= (2*pi*a*seglength);
    lamb = (1/wavelengthM)-(1/wavelength2M);
    NumPhotons = j*lamb*(1-g2);
  }
 
  cout<< NumPhotons<< endl;

  int iNumPhotons= (int)NumPhotons;
  cout << iNumPhotons<< endl;
  cout << endl;

  vector<double> xposition;
  vector<double> yposition;
 
  for (int i=0; i<20; i++){

    //get the wavelength
    double w = ran.Rndm();
    double wavelength_random = w*(lambda1inv-lambda2inv)+lambda2inv;
    double wavelengthPhoton =(1/wavelength_random);
  

    //get the photon altitude
    w=ran.Rndm();
    double photalt = w*(seglength)+(alt3M-seglength);
 

    //get the azimuthal angle
    w=ran.Rndm();
    double azim = w*360.0;
   

    //get our rotated segment


    TVector3 h(0,0,1);
    h.Print();
    cout << endl;

    TRotation r;
    r.RotateX(pi-theta);
    r.XX();
    r.YY();
    r.ZZ();
  
  
    TVector3 v=r*h;
    v.Print();
    cout << endl;
  
    TVector3 v1= v;
    TRotation r1;
    r1.RotateZ(azim*TMath::DegToRad());
  
    TVector3 v2=r1*v1;
  
    v2.Print();
    cout << endl;


    // get the direction of the photon
    TVector3 unitphoton;
    unitphoton.SetMagThetaPhi(1.,(pi-theta),azim);
    unitphoton.Print();
    TVector3 v3 = r*unitphoton;
    TVector3 v4 = r1*v3; 
    v4.Print();
    cout << endl;
    cout << endl;


    //get the distance the photon travels to the ground
    TVector3 k(0,0,1);
    double d = h.Dot(k)/v4.Dot(k);
    cout << d << endl;


    // get our position vector on the ground
    TVector3 positionvect = d*v4-h;
    positionvect.Print();

    xposition.push_back(positionvect.X());
    yposition.push_back(positionvect.Y());

  }

  cout << xposition.size() << endl;
  for (int i=0; i<xposition.size(); i++){
    
    cout << xposition[i] << "  " << yposition[i] << endl;
  }



      return 0; }






















 
