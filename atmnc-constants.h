#ifndef atmnc_constants_h__
#define atmnc_constants_h__

namespace Constants{
  double const rEarth=6378.140e3; // [m]
  // double const pi=3.14159265358979324; -> use TMath
  //  double const deg2rad=pi/180;        -> use TMath
  //  double const rad2deg=180./pi;       -> use TMath
  double const c=2.99792458e8;
  double const eps=1.e-6;
  double const Aair=14.5;
  double const mb2m2=1.e-31;
  double const abogN=6.0221367e26;
  double const rho2numDens=1./(Aair/abogN/mb2m2);
}


#endif // atmnc_constants_h__
