#ifndef ContainerClass_h__
#define ContainerClass_h__

#include <myheader/myRoot.h>
#include <myheader/myStandard.h>

int const MAX_ARRAY_SIZE = 32768;

struct ATMNCStatInfo{
  ATMNCStatInfo();
  void branch(TTree* tree);
  void setBranchAddress(TTree* tree);
  double crsamp;
  double crsim;
  double fintsum;
  //  double fintEach[4]; // element size is dummy. acturally it is [n1ry]
  double sarea;
  double omega;
  double tequiv;
  double fact1ry;
  static int const mtbl = 151; // hard-coded in include/atmnc-primary.inc
  double ebin1ry[2][mtbl]; //[before/after cutoff][n1ry][mtbl]
  double nobs1ry[2][4][mtbl]; // [4] is dummy, to avoid variable array size.
  double dflux1ry[2][4][mtbl];
};


struct ATMNCHistInfo{
  ATMNCHistInfo();
  void branch(TTree* tree);
  void setBranchAddress(TTree* tree);
  static int const nbin = 103; // hard-coded in include/
  int isite;
  int icosz;
  int iazim;
  double Darea;
  double Domega;
  double coszFrom;
  double coszTo;
  double azimFrom;
  double azimTo;

  int nkind;
  int kind[16];
  double ebin[nbin];
  double nobs[16][nbin];
  double nobsErr[16][nbin];
  double flux[16][nbin];
};



struct ATMNCInputInfo{
  ATMNCInputInfo();
  void init();
  void dump();
  void branch(TTree *tree);
  void setBranchAddress(TTree *tree);

  // for the simulation
  //  static int const fnameLength = 100;
  //  char pref[fnameLength];
  std::string* m_inputFile;
  std::string* pref;
  int n1ry;
  int kind1ry[16];          // [n1ry] kind1ry
  //  char f1ryFile[16][fnameLength]; // [n1ry] data table file
  std::string* f1ryFile[16];
  //  char cfil[fnameLength]; // cutoff table
  //  char gfil[fnameLength]; // geomagnetic field table (IGRF)
  //  char afil[fnameLength]; //  Air density table
  std::string* cfil; // cutoff table
  std::string* gfil; // geomagnetic field table (IGRF)
  std::string* afil; //  Air density table
  int mnth;   // start_month*100 + month_end
  double dmyr; // Year from IGRF date
  double epmi; // Min energy for primary CRs
  double lat0; // Center Lattitide of CR injection
  double lon0; // Center Longitude of CR injection
  double rsmp; // Sampling radius in degree
  double czmi; // Sampling minimum cos(zenith)
  double hinj; // CR injection Height
  double ehmi; // Minimum hadron int. energy
  double eemi; // Minimum elemag int. energy
  double eprp; // Minimum Energy for all particle propagation
  double tlim; // Maximum time for p,n,e,g propagation
  double rbnd; // Maximum distance from earth center
  double epht; // below this use jam int model
  bool noem;   // swith off elemag int.
  bool auti;  // Auto generation of random seed
  int nwup;  // Call rndc before simulation
  bool dbg1;  // Record random seed for each event
  bool dbg2;  // Record random seed for each particle
  int seed[2];  // Seeds in case for auti
  //char vfil[fnameLength]; //
  //  char vpfl[fnameLength]; //
  std::string* vfil;
  std::string* vpfl;
  int neve; // CR// to sample
  double fcrs; // Factor for cross section
  double rpft; // Production ratio for Jam
  double jobt; // Maximum Jobtime
  
  // for the observation
  //  char rfil[fnameLength];
  //  char ofil[fnameLength];
  std::string* rfil;
  std::string* ofil;
  bool is1d;
  int nsite;
  double siteParam[16][4]; //[nsite][nparam(=4)]
  double comi;  // Minimum cos(zenith) to obs
  int nczo;  // Observation cos(zenith) devision #
  int nazo;  // Observation azimuth devision #
  int nkind;
  int obsp[16]; //[nkind]
  double eomi; // Minimum energy to obs
  double etmi; // Minimum energy in table
  bool alpr; // Print when observed with 1ry
  int nstp[2]; // print every [0] events and stop at [1] prints
};


#endif // ContainerClass_h__
