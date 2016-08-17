//Wrapper for atmosphere class  as needed by kascade (fortran callable).
//Mainly matches types.

#include <cstdlib>
#include <iostream>

#include <string>
#include "atmosphere.h"
using namespace std;

#include "kasatmGrISU.h"

static atmosphere atm;

const int LAMBDA_MIN = 180;
const int LAMBDA_MAX = 900;
const int  STEP_SIZE_NEEDED = 5;
const int LSTEPS = ( 1 + (LAMBDA_MAX - LAMBDA_MIN) / STEP_SIZE_NEEDED ); 
const int ALT_STEPS = 51;
const float NM = 1.E-09;
const float C_LIGHT = (2.99792E+08F); // Speed of light in vacuum

void read_extint(char filename[],float extint[][LSTEPS]);
void run_error1(const char string[]);

void read_extint(char filename[],float extint[][LSTEPS])
{
  /* Program to read in extinction coefficients. Very finely tuned to file
     format in "kextint.dat"                                            */

  int    i,j,bin,mid_step,step,lambda,lambda_in;
  char   err_message[150];
  FILE   *finp;

  finp = fopen(filename, "r" );

  if( finp != NULL )
    {
      for(lambda = LAMBDA_MIN,step = STEP_SIZE_NEEDED ; 
          lambda <= LAMBDA_MAX; lambda += step )
        {
	  /* Read in wavelength. */
          fscanf(finp,"%d %*[^\n]\n",&lambda_in);

          if(lambda != lambda_in)
            run_error1("Error in reading kextint.dat");

          bin=(lambda - LAMBDA_MIN) / STEP_SIZE_NEEDED;

	  /* For given wavelength, read in  extinction coefficients.      */
          for(i=0 ;i<ALT_STEPS ;i++)
              fscanf(finp," %f %*c",&extint[i][bin]); 

	  /* F 
ill in blanks in coefficients array, if any */
          if(step != STEP_SIZE_NEEDED)
            {
              mid_step = step / STEP_SIZE_NEEDED;

              for(i=1 ;i<mid_step ;i++)
                for(j=0 ;j<ALT_STEPS ;j++)
                  extint[j][bin - mid_step + i] = extint[j][bin - mid_step] +
                    ( (float) i / mid_step ) * 
                      ( extint[j][bin] -  extint[j][bin - mid_step] );
            }

          switch( lambda )     /* These step sizes are specific to file 
                                  format in kextint.dat              */
            {
            case 270:
              step=10;
              break;

            case 280:
              step=20;
              break;

            case 400:
              step=50;
              break;

            case 700:
              step=100;
              break;
            }
        }
    }

  else
    {
      (void) sprintf( err_message, "Unable to open kextint data file: %80s ",
		filename );
      run_error1( err_message );
    }
  
  fclose(finp);
}
/***************** end of read_extint() ****************************/
void run_error1(const char string[])
{
  fprintf(stderr,"\n%s\n",string);
  exit(EXIT_FAILURE);
  
}
/***************** end of run_error1 ****************************/
float atm_prob(char filename[],float* zp, float* hobsp, 
               float* dnp, float* wavep)
{
  /* Returns 1 if the photon survives the atmosphere, zero if it does not.
     The probability of passing is determined by the optical depth, which is
     calculated from the altitude of emission, the observatory altitude and 
     the photon wavelength.                                               */

  int iwave,ihobs,ihgt;
  float tlow,thigh,atmprob;

  static float extint[ALT_STEPS][LSTEPS];
  static int read_flag = 0;
  float z,hobs,dn,wave;
  
  z = *zp;
  hobs = *hobsp;
  dn = *dnp;
  wave = *wavep;

  /* fill extint array if this is first call to atm_pass */
  if ( read_flag == 0) {
    read_flag = 1;
    read_extint(filename,extint);
  }

  /* Convert from m to nm. */

  wave /= NM;
  iwave = (int)( (wave - 180) / 5 );

  ihobs = (int)( hobs / 1000 );

  tlow = extint[ihobs][iwave] + ( extint[ihobs+1][iwave] - 
				 extint[ihobs][iwave] ) *
				   (hobs / 1000. - ihobs);

  ihgt = (int)(z / 1000);
  

  thigh = extint[ihgt][iwave] + ( extint[ihgt+1][iwave] -
				 extint[ihgt][iwave] ) *
				   (z / 1000. - ihgt);

  //cout << "wave,ihgt,ihobs,tlow,thigh " << wave << " " 
  //   << ihgt << " " << ihobs << " " << tlow << " " << thigh << endl;
  /* tlow and thigh are, respectively, the optical depths at the observation
     and emission altitudes.                                              */
  if(dn != 0.0)
    {
      atmprob = exp( (tlow - thigh) / dn );
     
    }
  else
      atmprob = 0.0;    

  return(atmprob);
}
/***************** end of atm_pass() *******************/

void initAtmosphere(string atmProfileFile)
// ***************************************************
// Initalize for the selected atmosphere.
// ***************************************************
{
  atm.init(atmProfileFile);
  return;
}
// ***************************************************

string getAtmosphereTitle()
// *************************************************************
//  Return the title of the atmosphere: "US76" or first line of atmprof file
// *************************************************************
{
  string atmTitle=atm.getTitle();
  return atmTitle;
}

// *************************************************************
// The following are wrappers for  calls from fortran, thus the arguments
// are pointers.
// **************************************************************

float rho(float* pAltitude)
  //altitude in meters,  density gm/cm**3
{
  float r = (float) atm.getDensity(pAltitude);
  return r;
}

float gms(float* pAltitude)
  //altitude in meters,  depth gm/cm**2
{
  float g =(float) atm.getDepth(pAltitude);
  return g;
}

 //altitude in meters,  depth gm/cm**2
float yds(float* pDepth)
{
  float y=(float) atm.getAltitudeM(pDepth);
  return y;
}

 //altitude in meters,  wavelength in nm
float eta(float* pAltitude, float* pLambdaNM)
{
  float e=(float) atm.getEta(pAltitude, pLambdaNM);
  return e;
}

float photon_time(float* e_heightp, float* d_heightp, 
                             float* zdircosp, float* wavep) 
{

/* procedure to calculate the time in nanoseconds for a photon of wavelength wave 
   to travel from an initial height e_height to a height d_height where the photon has a 
   z- direction cosine of zdircos.  The z-axis is pointed vertically upward.
   The procedure assumes that e_height > d_height and zdircos < 0.

*/

  float tsec, alpha;
  float delh;
  float delgms;

  float e_height,d_height,zdircos,wave;
  static float rho_sea;
  static int ifirst = 0;
  float zeroAlt = 0.0;

  if (!ifirst) {
    ifirst = 1;
    rho_sea = rho(&zeroAlt);
  }

  e_height = *e_heightp;
  d_height = *d_heightp;
  delh = fabs(e_height - d_height);
  zdircos = fabs(*zdircosp);
  wave = *wavep;

  wave = wave/NM;

  delgms = fabs(gms(e_heightp) - gms(d_heightp));
  
  alpha = eta(&zeroAlt,&wave)/rho_sea;

  float tnoair = delh/((C_LIGHT)*zdircos);
  
  tsec = alpha*delgms/( (C_LIGHT)*100.*zdircos);  // 100 gets units right

  tsec = tsec + tnoair;
  return(tsec/(NM));

}
/**************end of photon time*********************/
float index_of_ref(float* lambdap, float* heightp) {
  float wavelength = (*lambdap)/(NM);
  float altM = *heightp;
  float Eta = eta(&altM,&wavelength);
  float IndexRef = 1.0 + Eta;
  return IndexRef;

}
