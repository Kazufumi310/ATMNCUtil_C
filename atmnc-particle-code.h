#ifndef ParticleCode_h__
#define ParticleCode_h__

namespace ParticleCode {
  int const nParType = 22;
  enum{ kep=1, kem=2, kgamma=3,
	kkp=4, kkm=5, kk0l=6, kk0s=7,
	kneut=8, kneutbar=9, kpro=10, kprobar=11,
	kpip=12, kpim=13, kpi0=14,
	kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
	knumubar=20,
	kalpha=21 };     
  
  //     Particle masses
  double const aMk0=0.497672, aMkc=0.493677, aMpi0=0.1349764, 
    aMpic=0.13956995,aMmu=0.105658389, aMe=0.00051099907,
    aMnu=0., aMg=0., aMneut=0.93956563, aMpro=0.93827231,
    aMalpha=3.72737;
  
  double const mass[nParticle]={ 
    //     *  e-   e+  gam, 
    0, aMe, aMe, aMg, 
    //     K+,   K-    KL    Ks
    aMkc, aMkc, aMk0, aMk0,
    //     N       N-bar   Proton Proton-bar 
    aMneut, aMneut, aMpro, aMpro,
    //     pi+    pi-    pi0
    aMpic, aMpic, aMpi0,    
    //     mu+   mu-   nue   a-nue numu  a-numu
    aMmu, aMmu, aMnu, aMnu, aMnu, aMnu, 
    //     alpha
    aMalpha };
}

#endif // ParticleCode_h__
