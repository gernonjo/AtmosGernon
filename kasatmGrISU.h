
/*
  New atmospheric code that uses the corsika atmospheric files,
  included the 21 and 22 atmospheres for Veritas.

  All atmospheric functions including photon_time and photon extinction
  are included here. The extinction code uses the usual kextint.dat and 
  the V.Summer.US76.50km.profileGrISU.ext and 
  V.Winter.US76.50km.profileGrISU.ext files for Veritas.  These extinction
  files are extended to 180 through 900 nanometers using the low and
  high parts of kextinct.dat. Limiting the wavelength range in cherenkf7 to 
  200 through 900 NM thus ignores the low and high table entries of 
  kextinct.dat. In so doing no code changes are required in cherenkf7; yet
  the atmosphere now agree with those used in the Corsika and in Glenn's
  KASCADE simulations.

 */

extern "C" {
  // returns atmospheric density in gm/cm3
  // pAltitude in meters
  float rho(float* pAltitude);//Use of pointers so that these are

  // returns thickness in gm/cm2
  // pAltitude in meters
  float gms(float* pAltitude);// fortran callable

  // returns altitude in meters
  // thickness in gm/cm2
  float yds(float* pDepth);

  // return 1 - index of refraction
  // pAltitude in meters,  pLambdaNM in nanometers
  float eta(float* pAltitude, float* pLambdaNM);
  
  // lambdap in METERS, heightp in meters
  float index_of_ref(float* lambdap, float* heightp);
  
  // initialize atmosphere
  // atmProfileFile example: atmprof6.dat (Corsika USStd Atmosphere)
  void initAtmosphere(string atmProfileFile);

  // returns title of atmosphere. First line of atmprofxxx.dat file
  string getAtmosphereTitle();
  
  // returns probability of photon passing through the atmosphere
  // filename: extinction file, read during the first call to this function
  // zp emission height in meters, hobsp observatory height in meters
  // dnp z direction cosine, dnp > 0.0. (downward pointing z axis)
  // wavep in METERS
  float atm_prob(char filename[],float* zp, float* hobsp, 
                 float* dnp, float* wavep);
  

  // returns photon time in nsec,
  // e_heightp emission height in meters
  // d_heightp detection height in meters
  // zdircosp z direction cosine > 0
  // wavep in METERS
  float photon_time(float* e_heightp, float* d_heightp, 
                             float* zdircosp, float* wavep);
}
